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
import itertools

from urllib.parse import urlparse, urlunparse

from playhouse import db_url as playhouse_db_url

from . import schema

_sqlite_adapter = None


def get_default_sqlite_adapter(quiet, apsw_warning=True):
    global _sqlite_adapter
    if _sqlite_adapter is None:
        _sqlite_adapter = _get_default_sqlite_adapter(quiet, apsw_warning=apsw_warning)
    return _sqlite_adapter


def _get_default_sqlite_adapter(quiet, apsw_warning):
    # Before mmpdb 3.0 we recommended using apsw because that allowed
    # copying the database into memory, which resulted in faster
    # searches. Python 3.7 added sqlite.Connection.backup, which
    # we can now assume exists for all users.
    return "sqlite"

def open_as_schema_database(playhouse_db):
    c = playhouse_db.cursor()
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
    if mmpdb_version != 4:
        raise DBError("Expecting mmpdb version 4, not %d" % (mmpdb_version,))

    import peewee
    if isinstance(playhouse_db, peewee.PostgresqlDatabase):
        return schema.PostgresMMPDatabase(playhouse_db)
    else:
        return schema.MMPDatabase(playhouse_db)


class DBError(Exception):
    def __init__(self, text):
        self.text = text

    def __str__(self):
        return self.text

    def __repr__(self):
        return "DBError(%r)" % (self.text,)


def _sqlite_copy_to_memory(db, quiet):
    import sqlite3
    
    db_connection = db.connection()
    if not isinstance(db_connection, sqlite3.Connection):
        sys.stderr.write(
            "WARNING: 'copy_to_memory' requires a SQLite database using Python's built-in sqlite connection, "
            f"not {db.__class__.__name__}. "
            "Keeping on-disk.\n"
            )
        return db

    memory_db = db.__class__(":memory:")
    
    import time
    
    if not quiet:
        t1 = time.time()
        sys.stderr.write("Copying database to memory...")
        sys.stderr.flush()

    # available starting in Python 3.7
    db_connection.backup(memory_db.connection())
        
    if not quiet:
        t2 = time.time()
        sys.stderr.write("\rDatabase copy took %.1f seconds.\n" % (t2 - t1,))

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

    def open_database(self, *, copy_to_memory=False, quiet=False, apsw_warning=True):
        if not os.path.exists(self.name):
            raise DBError("File does not exist")
        database_class = playhouse_db_url.schemes[
            get_default_sqlite_adapter(quiet, apsw_warning=apsw_warning)
            ]
        try:
            db = database_class(database=self.name)
        except Exception as err:
            raise DBError(str(err))
        if copy_to_memory:
            db = _sqlite_copy_to_memory(db, quiet)
        return open_as_schema_database(db)


class DBUrl(DBInfo):
    def __repr__(self):
        return "DBUrl(%r)" % (self.name,)

    def get_human_name(self):
        return "database url %r" % (self.name,)

    def open_database(self, *, copy_to_memory=False, quiet=False, apsw_warning=True):
        try:
            db = playhouse_db_url.connect(self.name)
        except Exception as err:
            raise DBError(str(err))
        if copy_to_memory:
            db = _sqlite_copy_to_memory(db, quiet)
        return open_as_schema_database(db)


def get_database_server(url):
    parsed = urlparse(url)
    if parsed.scheme not in playhouse_db_url.schemes:
        return None

    if parsed.scheme == "postgres":
        # Are we pointing to the server, or a specific database
        if parsed.path in ("", "/", "/postgres"):
            new_url = urlunparse(parsed._replace(path="/postgres"))
            db = playhouse_db_url.connect(new_url)
            return PostgresServer(url, db)
        return None

    # I don't know how to support other database servers
    return None
    
class PostgresServer(object):
    def __init__(self, url, db):
        self.url = url
        self.db = db
        self._parsed = urlparse(self.url)
        
    def get_mmpdb_databases(self, reporter):
        parsed = self._parsed
        c = self.db.cursor()
        
        # List the available databases.
        c.execute("SELECT datname FROM pg_database WHERE datistemplate = false")
        databases = [row[0] for row in c if row[0] != "postgres"]
        
        # See which has the right schema
        for database in databases:
            if self.is_mmpdb_database(database, reporter):
                url = urlunparse(parsed._replace(path="/" + database))
                yield DBUrl(url)

    def is_mmpdb_database(self, database, reporter):
        import psycopg2.errors
        url = urlunparse(self._parsed._replace(path = "/" + database))
        # XXX need error checking
        db = playhouse_db_url.connect(url)
        db.connect()
        c = db.cursor()
        try:
            c.execute("SELECT mmpdb_version FROM dataset LIMIT 2")
        except psycopg2.errors.UndefinedTable:
            return False
        values = [row[0] for row in c]
        if not values:
            return True

        if len(values) > 1:
            # Shouldn't happen
            return False

        mmpdb_version = values[0]
        if mmpdb_version == 4:
            return True
        reporter.warning(f"Skipping unsupported version {mmpdb_version} database database: {url!r}")
        return False
        
    
def is_valid_dburl(url):
    parsed = urlparse(url)
    if (parsed.scheme in playhouse_db_url.schemes):
        return True
    if parsed.scheme == "apsw":
        sys.stderr.write("Cannot use `apsw:` databases - apsw module not available.")

    return False


def get_dbinfo(dburl):
    # If you have a filename which looks like a database url
    # then be explicit and use "apsw:" or "sqlite:".
    # (XXX add new scheme to autodetect the best supported API?)
    if is_valid_dburl(dburl):
        return DBUrl(dburl)
    return DBFile(dburl)


# just the filenames. Does not join with the dirname
def _get_mmpdb_filenames(dirname):
    return sorted(filename for filename in os.listdir(dirname) if filename.lower().endswith(".mmpdb"))


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
        server = get_database_server(database)
        if server is not None:
            for entry in server.get_mmpdb_databases(reporter):
                yield entry
            return
                
        if os.path.isdir(database):
            for filename in get_mmpdb_filenames_in_directory(database):
                yield DBFile(filename)

        elif is_valid_dburl(database):
            yield DBUrl(database)

        elif os.path.isfile(database):
            yield DBFile(database)

        else:
            reporter.report(f"Not a file, directory, or supported database URL:{database!r}")


def iter_dbinfo_and_dataset(databases, reporter, apsw_warning=True):
    for dbinfo in iter_dbinfo(databases, reporter):
        reporter.update("Opening %s ... " % (dbinfo.get_human_name(),))
        database = None
        try:
            database = dbinfo.open_database(
                quiet=reporter.quiet,
                apsw_warning=apsw_warning,
                )
            dataset = database.get_dataset()
        except DBError as err:
            reporter.update("")
            reporter.report("Skipping %s: %s" % (dbinfo.get_human_name(), err))
            if database is not None:
                database.close()
            continue
        reporter.update("")
        try:
            yield dbinfo, dataset
        finally:
            database.close()


def open_database(dburl, copy_to_memory=False, quiet=False, apsw_warning=True):
    return get_dbinfo(dburl).open_database(
        copy_to_memory=copy_to_memory,
        quiet=quiet,
        apsw_warning=apsw_warning,
        )


#####


def reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id, cursor, reporter):
    # Mapping from rule environment id to rule environment statistics id

    reporter.update("Computing aggregate statistics")
    num_pairs = dataset.get_num_pairs(cursor=cursor)

    all_pairs = dataset.iter_pairs(cursor=cursor)
    all_pairs = reporter.progress(all_pairs, "Computing and updating aggregate statistics", num_pairs)

    def generate_stats():
        for rule_environment_id, rule_environment_pairs in itertools.groupby(
            all_pairs, (lambda pair: pair.rule_environment_id)
        ):

            rule_environment_pairs = list(rule_environment_pairs)  # now a list, not iterator

            for property_name_id in property_name_ids:
                deltas = []
                compound_values = compound_values_for_property_name_id[property_name_id]
                for pair in rule_environment_pairs:
                    value1 = compound_values.get(pair.compound1_id, None)
                    if value1 is None:
                        continue
                    value2 = compound_values.get(pair.compound2_id, None)
                    if value2 is None:
                        continue
                    deltas.append(value2 - value1)
                if deltas:
                    from . import index_algorithm

                    stats = index_algorithm.compute_aggregate_values(deltas)
                    yield (rule_environment_id, property_name_id, stats)

    stats_info = generate_stats()

    reporter.update("Getting information about which rule statistics exist...")
    existing_stats_ids = dataset.get_rule_environment_statistics_mapping(property_name_ids, cursor=cursor)

    seen_stats_ids = set()
    num_updated = num_added = 0
    for (rule_environment_id, property_name_id, stats) in stats_info:
        key = (rule_environment_id, property_name_id)
        stats_id = existing_stats_ids.get(key, None)
        if stats_id is not None:
            dataset.update_rule_environment_statistics(stats_id, stats)
            seen_stats_ids.add(stats_id)
            num_updated += 1
        else:
            dataset.add_rule_environment_statistics(rule_environment_id, property_name_id, stats)
            num_added += 1

    to_delete = set(existing_stats_ids.values()) - seen_stats_ids
    num_deleted = len(to_delete)
    if to_delete:
        delete_report = reporter.progress(iter(to_delete), "Deleting rule statistics", num_deleted)
        while 1:
            ids = list(itertools.islice(delete_report, 0, 1000))
            if not ids:
                break
            dataset.delete_rule_environment_statistics(ids)

    reporter.report(
        "Number of rule statistics added: %d updated: %d deleted: %d" % (num_added, num_updated, num_deleted)
    )
