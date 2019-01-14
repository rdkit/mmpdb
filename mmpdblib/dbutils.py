# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

# Commands to work with database urls and figure out which are local
# database files.

## NOTE: Much of this comes from an older version of mmpdb which
## supported MySQL *and* SQLite. This version of mmpdb only supports
## SQLite. I have not yet cleaned up the dead code.

import os
import sys

try:
    from urlparse import urlparse  # Python 2
except ImportError:
    from urllib.parse import urlparse  # Python 3

from .playhouse import db_url as playhouse_db_url

from . import schema

_sqlite_adapter = None
def get_default_sqlite_adapter(quiet):
    global _sqlite_adapter
    if _sqlite_adapter is None:
        _sqlite_adapter = _get_default_sqlite_adapter(quiet)
    return _sqlite_adapter


def _get_default_sqlite_adapter(quiet):
    try:
        import apsw
        return "apsw"
    except ImportError:
        if not quiet:
            sys.stderr.write("WARNING: APSW not installed. Falling back to Python's sqlite3 module.\n")
        return "sqlite"
        

def open_as_schema_database(playhouse_db):
    c = playhouse_db.get_cursor()
    try:
        c.execute("SELECT id, mmpdb_version from dataset LIMIT 2")
    except Exception:
        raise DBError("Missing required 'dataset' table")
    values = list(c)
    if not values:
        raise DBError("Does not contain a dataset")
    if len(values) > 1:
        raise DBError("Contains too many datasets")
    id, mmpdb_version = values[0]
    if id != 1:
        raise DBError("The dataset has the wrong id")
    if mmpdb_version != 2:
        raise DBError("Expecting mmpdb version 2, not %d" % (mmpdb_version,))
    
    return schema.MMPDatabase(playhouse_db)
    

class DBError(Exception):
    def __init__(self, text):
        self.text = text

    def __str__(self):
        return self.text

    def __repr__(self):
        return "DBError(%r)" % (self.text,)


def _apsw_copy_to_memory(db, quiet):
    try:
        import apsw
    except ImportError:
        sys.stderr.write("WARNING: 'copy_to_memory' requires the apsw module. Keeping on-disk.\n")
        return db

    from .playhouse.apsw_ext import APSWDatabase
    if not isinstance(db, APSWDatabase):
        sys.stderr.write("WARNING: 'copy_to_memory' requires an APSW database, not '%s'. Keeping on-disk.\n"
                         % (db.__class__,))
        return db

    disk_conn = db.get_conn()  # the low-level APSW connection
    
    memory_db = APSWDatabase(":memory:", **db.connect_kwargs)
    memory_conn = memory_db.get_conn()  # the low-level APSW connection

    import time
    if not quiet:
        t1 = time.time()
        sys.stderr.write("Copying database to memory...")
        sys.stderr.flush()
    with memory_conn.backup("main", disk_conn, "main") as backup:
        backup.step()
    if not quiet:
        t2 = time.time()
        sys.stderr.write("\rDatabase copy took %.1f seconds.\n" % (t2-t1,))

    return memory_db

    
# Base class
class DBInfo(object):
    def __init__(self, name):
        self.name = name

    def get_human_name(self):
        raise NotImplementedError

    def open_database(self):
        raise NotImplementedError


class DBFile(DBInfo):
    def __repr__(self):
        return "DBFile(%r)" % (self.name,)

    def get_human_name(self):
        return "file %r" % (self.name,)

    def open_database(self, copy_to_memory=False, quiet=False):
        if not os.path.exists(self.name):
            raise DBError("File does not exist")
        database_class = playhouse_db_url.schemes[get_default_sqlite_adapter(quiet)]
        try:
            db = database_class(database=self.name)
        except Exception as err:
            raise DBError(str(err))
        if copy_to_memory:
            db = _apsw_copy_to_memory(db, quiet)
        return open_as_schema_database(db)
    

class DBUrl(DBInfo):
    def __repr__(self):
        return "DBUrl(%r)" % (self.name,)

    def get_human_name(self):
        return "database url %r" % (self.name,)
    
    def open_database(self, copy_to_memory, quiet=False):
        try:
            db = playhouse_db_url.connect(self.name)
        except Exception as err:
            raise DBError(str(err))
        if copy_to_memory:
            db = _apsw_copy_to_memory(db)
        return open_as_schema_database(db)


def is_valid_dburl(url):
    parsed = urlparse(url)
    return parsed.scheme in playhouse_db_url.schemes


def get_dbinfo(dburl):
    # If you have a filename which looks like a database url
    # then be explicit and use "apsw:" or "sqlite:".
    # (XXX add new scheme to autodetect the best supported API?)
    if is_valid_dburl(dburl):
        return DBUrl(dburl)
    return DBFile(dburl)
    

# just the filenames. Does not join with the dirname
def _get_mmpdb_filenames(dirname):
    return [filename for filename in os.listdir(dirname) if filename.lower().endswith(".mmpdb")]    

# files in ".", without the leaing "./"
def get_mmpdb_filenames_in_current_directory():
    return _get_mmpdb_filenames(".")

# Files including the path to the directory
def get_mmpdb_filenames_in_directory(dirname):
    return [os.path.join(dirname, filename) for filename in _get_mmpdb_filenames(dirname)]

# Get a list of all of the databases:
#   - if it's a directory, look for *.mmpdb files in that directory
#   - if it looks like a databse URL then use playhouse to open it
#   - if it's a file, use apsw/sqlite to open it
#   - otherwise, no clue.

def iter_dbinfo(databases, reporter):
    if not databases:
        for filename in get_mmpdb_filenames_in_current_directory():
            yield DBFile(filename)
        return

    for database in databases:
        if os.path.isdir(database):
            for filename in get_mmpdb_filenames_in_directory(database):
                yield DBFile(filename)

        elif is_valid_dburl(database):
            yield DBUrl(database)

        elif os.path.isfile(database):
            yield DBFile(database)

        else:
            reporter.report("Not a file, directory, or supported database URL: %r"
                            % (database,))


def open_database(dburl, copy_to_memory=False, quiet=False):
    return get_dbinfo(dburl).open_database(copy_to_memory=copy_to_memory, quiet=quiet)


def open_database_from_args_or_exit(args):
    dburl = args.databases[0]
    dbinfo = get_dbinfo(dburl)
    copy_to_memory = getattr(args, "in_memory", False)
    quiet = getattr(args, "quiet", False)
    try:
        return dbinfo.open_database(copy_to_memory=copy_to_memory, quiet=quiet)
    except DBError as err:
        sys.stderr.write("Cannot connect to %s: %s\n"
                         % (dbinfo.get_human_name(), err))
        raise SystemExit(1)


def open_dataset_from_args_or_exit(args):
    db = open_database_from_args_or_exit(args)
    return db.get_dataset()
