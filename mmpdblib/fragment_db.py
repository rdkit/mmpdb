# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2021, Andrew Dalke Scientific, AB
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

import os
import sqlite3
import dataclasses
import importlib.resources

from .fragment_types import (
    FragmentOptions,
    FragmentRecord,
    FragmentErrorRecord,
    Fragmentation,
)
from . import schema

#### "fragdb" -- SQLite-based fragment file

# NOTE: There is configuration information in three files!
# 1) fragment_types.py -- the data types
# 2) fragment_schema.sql -- defines the SQL schema
# 3) fragment_db.py -- (this file) defines the mapping from SQL to the data types


_schema_template = None


def get_schema_template():
    global _schema_template
    if _schema_template is None:
        _schema_template = importlib.resources.read_text(
            "mmpdblib",
            "fragment_schema.sql",
        )
    return _schema_template

def get_fragment_create_index_sql():
    return importlib.resources.read_text(
            "mmpdblib",
            "fragment_create_index.sql",
        )

def init_fragdb(c, options):
    # Ensure SQL can read the file, and that no records exist
    try:
        c.execute("SELECT count(*) FROM record")
    except sqlite3.OperationalError as err:
        if "no such table: record" in str(err):
            pass
        else:
            raise
    else:
        raise AssertionError("I tried to delete the fragdb file but apparently it's valid?!")

    # To try:
    #  c.execute("PRAGMA journal_mode=WAL")
    c.execute("PRAGMA synchronous=off")
    # My WAL attempt left a couple of temp files in the local directory.
    # Perhaps I didn't close things?

    # Create the schema
    schema._execute_sql(c, get_schema_template())

    insert_options(c, options)


def open_fragdb(filename):
    with open(filename, "rb"):  # Let Python raise any IOErrors
        pass
    db = sqlite3.connect(filename)
    c = db.cursor()
    # See if this is a database
    try:
        c.execute("SELECT COUNT(*) FROM record")
    except sqlite3.OperationalError:
        raise ValueError(f"{filename!r} does not appear to contain a fragdb database (missing 'record' table?)")
    try:
        c.execute("SELECT COUNT(*) FROM options")
    except sqlite3.OperationalError:
        raise ValueError(f"{filename!r} does not appear to contain a fragdb database (missing 'options' table?)")

    for (count,) in c:
        if count != 1:
            raise ValueError(f"{filename!r} does not appear to contain a fragdb database (no options?)")
        break

    try:
        options = select_options(c)
    except sqlite3.OperationalError:
        raise ValueError(f"{filename!r} does not appear to contain a fragdb database (cannot read options?)")
    return FragDB(None, options, db, c)


#### options

_option_attrs = [field.name for field in dataclasses.fields(FragmentOptions)]
_option_attrs_v3 = [s for s in _option_attrs if s != "min_heavies_total_const_frag"]
_option_attrs_v4 = _option_attrs

_option_cols = ["version"] + _option_attrs

_option_fields = ", ".join(_option_cols)
_option_qs = ",".join("?" * len(_option_cols))

_insert_options_sql = f"INSERT INTO options({_option_fields}) VALUES({_option_qs})"


def insert_options(c, options):
    values = [4] + [getattr(options, attr) for attr in _option_attrs]
    c.execute(_insert_options_sql, values)


def select_options(c):
    for (version,) in c.execute("SELECT version FROM options"):
        if version == 3:
            attrs = _option_attrs_v3
        elif version == 4:
            attrs = _option_attrs_v4
        else:
            raise ValueError("Expected a version 3 or 4 database, not version {version}")
        break
    else:
        raise ValueError("Unable to determine the version")

    option_fields = ", ".join(attrs)
    sql = f"SELECT {option_fields} FROM options"
    for values in c.execute(sql):
        kwargs = dict(zip(attrs, values))
        return FragmentOptions(**kwargs)
    raise AssertionError("Missing options in fragdb")


#### FragmentRecord

_record_attrs = ("id", "input_smiles", "num_normalized_heavies", "normalized_smiles")
_record_cols = ("title", "input_smiles", "num_normalized_heavies", "normalized_smiles")
_record_fields = ", ".join(_record_cols)
_record_qs = ",".join("?" * len(_record_cols))

_insert_record_sql = f"INSERT INTO record (id, {_record_fields}) VALUES (?, {_record_qs})"


def insert_fragment_record(c, rec, record_id):
    record_values = [getattr(rec, attr) for attr in _record_attrs]
    c.execute(_insert_record_sql, [record_id] + record_values)

    for frag in rec.fragmentations:
        c.execute(
            _insert_fragmentation_sql,
            [record_id] + [getattr(frag, attr) for attr in _fragmentation_attrs],
        )


_select_record_by_title_sql = f"SELECT id,{_record_fields} FROM record WHERE title = ?"


def select_fragment_record_by_title(c, title):
    c.execute(_select_record_by_title_sql, (title,))
    for row in c:
        break
    else:
        return None

    record_id, *record_values = row

    return FragmentRecord(
        *record_values,
        fragmentations=list(select_fragmentations_by_record_id(c, record_id)),
    )


_fragmentation_attrs = [field.name for field in dataclasses.fields(Fragmentation)]
_fragmentation_cols = _fragmentation_attrs[:]
_fragmentation_fields = ",".join(_fragmentation_cols)
_fragmentation_qs = ",".join("?" * len(_fragmentation_cols))

_insert_fragmentation_sql = f"""
INSERT INTO fragmentation (record_id, {_fragmentation_fields}) VALUES (?, {_fragmentation_qs})
"""

_select_fragmentation_by_record_id_sql = f"SELECT {_fragmentation_fields} FROM fragmentation WHERE record_id = ?"


def select_fragmentations_by_record_id(c, record_id):
    c.execute(_select_fragmentation_by_record_id_sql, (record_id,))
    for row in c:
        yield Fragmentation(*row)


_select_fragmentations_sql = f"SELECT id,{_record_fields} FROM record"


def iter_fragment_records(record_c, fragmentation_c):
    record_c.execute(_select_fragmentations_sql)
    for row in record_c:
        record_id, *record_values = row
        fragmentations = list(select_fragmentations_by_record_id(fragmentation_c, record_id))
        yield FragmentRecord(
            *record_values,
            fragmentations=fragmentations,
        )


_select_fragmentations_error_sql = "SELECT title, input_smiles, errmsg FROM error_record"


def iter_fragment_error_records(record_c, fragmentation_c):
    record_c.execute(_select_fragmentations_error_sql)
    for row in record_c:
        title, input_smiles, errmsg = row
        yield FragmentErrorRecord(
            id=title,
            input_smiles=input_smiles,
            errmsg=errmsg,
        )


#### FragmentErrorRecord

_error_record_attrs = ("id", "input_smiles", "errmsg")
_error_record_cols = ("title", "input_smiles", "errmsg")
_error_record_fields = ", ".join(_error_record_cols)
_error_record_qs = ",".join("?" * len(_error_record_cols))

_insert_fragment_error_sql = f"""
INSERT INTO error_record ({_error_record_fields}) VALUES ({_error_record_qs})
"""


def insert_fragment_error_record(c, rec):
    c.execute(_insert_fragment_error_sql, [getattr(rec, attr) for attr in _error_record_attrs])


_select_error_record_by_title_sql = f"SELECT {_error_record_fields} FROM error_record WHERE title = ?"


def select_fragment_error_record_by_title(c, title):
    c.execute(_select_error_record_by_title_sql, (title,))
    for row in c:
        return FragmentErrorRecord(*row)

    return None


####


class FragDBWriter:
    def __init__(self, filename, db, c, options):
        self.filename = filename
        self.db = db
        self.c = c
        self._record_id = 0  # Manage the record ids myself

    def close(self):
        db, c = self.db, self.c
        if db is not None:
            schema._execute_sql(c, get_fragment_create_index_sql())
            self.db = self.c = None
            c.close()
            db.commit()  # don't rollback - keep partial writes too.
            db.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def cursor(self):
        return self.db.cursor()

    def write_version(self):
        raise NotImplementedError("cannot write version twice")

    def write_options(self, options):
        raise NotImplementedError("cannot add options twice")

    def write_records(self, fragment_records):
        c = self.c
        for rec in fragment_records:
            if rec.errmsg is not None:
                insert_fragment_error_record(c, rec)
                continue

            self._record_id = record_id = self._record_id + 1
            insert_fragment_record(c, rec, record_id)


###


class FragDB:
    def __init__(self, metadata, options, db, c):
        self.metadata = metadata
        self.db = db
        self.c = c
        self.options = options

    def get(self, id):
        obj = select_fragment_record_by_title(self.c, id)
        if obj is not None:
            return obj

        return select_fragment_error_record_by_title(self.c, id)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        c, db = self.c, self.db
        if self.c is not None:
            self.c = self.db = None
            c.close()
            db.rollback()

    def __iter__(self):
        return iter_fragment_records(self.db.cursor(), self.db.cursor())

    def iter_error_records(self):
        return iter_fragment_error_records(self.db.cursor(), self.db.cursor())

    def cursor(self):
        return self.db.cursor()

### connect to the database and prepare for writing.


def open_fragment_writer(filename, options):
    assert filename is not None
    # Remove any existing file.
    try:
        os.unlink(filename)
    except FileNotFoundError:
        pass
    db = sqlite3.connect(filename)
    c = db.cursor()
    writer = FragDBWriter(filename, db, c, options)
    init_fragdb(c, options)
    return writer
