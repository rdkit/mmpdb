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

from __future__ import print_function, absolute_import

import sys
import json
import re
import itertools
import os
import sqlite3
import dataclasses

from . import __version__ as mmpdblib_version

from . import config
from . import fileio
from . import fragment_algorithm
from . import reporters
from .fragment_types import (FragmentOptions, FragmentRecord, FragmentErrorRecord,
                                 Fragmentation, FragmentFormatError)
from ._compat import basestring
from . import schema

SOFTWARE = "mmpdb-" + mmpdblib_version

##### Read

parse_max_heavies_value = config.positive_int_or_none
parse_max_rotatable_bonds_value = config.positive_int_or_none
parse_min_heavies_per_const_frag = config.nonnegative_int


def parse_num_cuts_value(value):
    if value not in ("1", "2", "3"):
        raise ValueError("must be '1', '2', or '3'")
    return int(value)


def parse_method_value(value):
    if value not in ("chiral",):
        if value in ("hussain", "dalke"):
            raise ValueError("'chiral' is supported in mmpdb v2, not %r" % (value,))
        raise ValueError("must be 'chiral'")
    return value


class FragmentReader(object):
    def __init__(self, metadata, options, reader, location):
        self.version = metadata["version"]
        self.software = metadata["software"]
        
        self.options = options
        self._reader = reader
        self.location = location
        
    def __iter__(self):
        if self._reader is None:
            raise ValueError("I/O operation on closed file")
        return self._reader
    
    def __next__(self):
        if self._reader is None:
            raise ValueError("I/O operation on closed file")
        return next(self._reader)
    
    next = __next__

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
    
    def close(self):
        reader = self._reader
        if reader is None:
            return
        self._reader = None
        reader.close()
    
    
def read_fragment_records(source):
    assert source is not None
    lc = source.lower()
    if not (lc.endswith("fragments") or lc.endswith("fragments.gz")):
        return open_fragdb(source)

    if source is None:
        infile = fileio.open_input(None)
        filename = "<stdin>"
        close = None
        
    elif isinstance(source, basestring):
        infile = fileio.open_input(source)
        filename = source
        close = infile.close

    else:
        infile = source
        filename = getattr(source, "name", "<unknown>")
        close = None

    location = fileio.Location(filename)
    location.save(record_format="fragment")

    line_reader = enumerate(infile, 1)
    line_reader, metadata, options, options_dict = _get_options(line_reader, location)

    reader = _read_fragment_records(line_reader, close, location, options_dict)
    x = next(reader)
    assert x == "ready"

    return FragmentReader(metadata, options, reader, location)


_json_loads = None
_json_module_name = None


def get_json_loads():
    global _json_loads, _json_module_name
    if _json_loads is None:
        # Timings reported for a fragment file with 37177 lines
        # (35634 "RECORD" and 1534 "IGNORE" records.)
        try:
            # 40.05 seconds
            import ujson
            _json_loads = ujson.decode
            _json_module_name = "ujson"
        except ImportError:
            try:
                # 41.85 seconds
                import cjson
                _json_loads = cjson.decode
                _json_module_name = "cjson"
            except ImportError:
                # 55.5 seconds
                _json_loads = json.loads
                _json_module_name = "json"

    return _json_loads


_option_parser = {
    "cut_smarts": str,
    "max_heavies": parse_max_heavies_value,
    "max_rotatable_bonds": parse_max_rotatable_bonds_value,
    "method": parse_method_value,
    "num_cuts": parse_num_cuts_value,
    "rotatable_smarts": str,
    "salt_remover": str,
    "min_heavies_per_const_frag": parse_min_heavies_per_const_frag
    }


def _get_options(line_reader, location):
    options_dict = {}
    options = FragmentOptions(**config.DEFAULT_FRAGMENT_OPTIONS.to_dict())
    version = None
    software = None
    lineno = 0
    loads = get_json_loads()
    for lineno, line in line_reader:
        try:
            fields = loads(line)
        except ValueError as err:
            if lineno == 1:
                raise FragmentFormatError("The input does not appear to be a fragment file",
                                          location)
            raise
        #fields = line.rstrip("\n").split("\t")
        if version is None:
            if len(fields) == 2 and fields[0] == "VERSION":
                version = fields[1]
                if version != "mmpdb-fragment/2":
                    location.save(lineno=lineno)
                    raise FragmentFormatError("This reader only supports version 'mmpdb-fragment/2', not version %r"
                                              % (version,), location)
            else:
                location.save(lineno=lineno)
                raise FragmentFormatError("Missing VERSION from first line in the file", location)
            continue

        if len(fields) == 2 and fields[0] == "SOFTWARE":
            software = fields[1]
            continue
                    
        if len(fields) != 3 or fields[0] != "OPTION":
            # Push back the current line
            line_reader = itertools.chain([(lineno, line)], line_reader)
            return line_reader, {"version": version, "software": software}, options, options_dict

        _, name, value_str = fields
        if name not in _option_parser:
            location.save(lineno=lineno)
            raise FragmentFormatError("Unknown OPTION %r" % (name,), location)

        if name in options_dict and options_dict[name][1] != value_str:
            location.save(lineno=lineno)
            raise FragmentFormatError("OPTION %s already set to %r on line %d"
                                      % (name, options_dict[name][1], options_dict[name][0]),
                                      location)
        
        parser = _option_parser[name]
        try:
            value = parser(value_str)
        except ValueError as err:
            location.save(lineno=lineno)
            raise FragmentFormatError("Cannot understand option %s (%r): %s"
                                       % (name, value_str, err), location)

        setattr(options, name, value)
        options_dict[name] = (lineno, value_str)
        
    # No remaining data
    if version is None:
        location.save(lineno=lineno)
        raise FragmentFormatError("Missing required VERSION line")
    
    return iter([]), {"version": "version", "software": software}, options, options_dict
            

def _read_fragment_records(line_reader, close, location, options_dict):
    recno = 0
    lineno = 0
    line = None

    def get_recno():
        return recno

    def get_lineno():
        return lineno

    def get_record():
        return line

    location.register(get_recno=get_recno,
                      get_lineno=get_lineno,
                      get_record=get_record,
                      )
    
    loads = get_json_loads()
    Fragmentation = fragment_algorithm.Fragmentation
    yield "ready"
    try:
        for lineno, line in line_reader:
            try:
                fields = loads(line)
            except ValueError as err:
                err.message = err.message.replace("line 1", "line %d" % (lineno,))
                err.args = (err.message,) + err.args[1:]
                raise
            label = fields[0]
            if label == "RECORD":
                assert label == "RECORD"
                recno += 1
                try:
                    _, id, input_smiles, num_normalized_heavies, normalized_smiles, fragment_fields_list = fields
                except ValueError:
                    raise FragmentFormatError("Expected 7 fields on RECORD line, not %d" % (len(fields),),
                                                   location)

                try:
                    # This is the hot spot for the reader. About 70% of the time is spent here.
                    fragmentations = [Fragmentation(*fragment_fields) for fragment_fields in fragment_fields_list]
                except TypeError:
                    # Try to report the error a bit more nicely:
                    for fragment_i, fragment_fields in enumerate(fragment_fields_list):
                        if len(fragment_fields) != 10:
                            raise FragmentFormatError("Expected fragment[%d] with 10 fields, not %d (%r)"
                                                           % (fragment_i, len(fragment_fields), fragment_fields),
                                                           location)
                    raise AssertionError  # I don't know what caused this error

                yield FragmentRecord(
                    id, input_smiles, num_normalized_heavies, normalized_smiles, fragmentations)
                continue

            if label == "IGNORE":
                try:
                    _, id, input_smiles, errmsg = fields
                except ValueError:
                    raise FragmentFormatError("Expected 4 fields on IGNORE line, not %d" % (len(fields),),
                                              location)
                yield FragmentErrorRecord(id, input_smiles, errmsg)
                continue
            
            ### It looks like this code isn't used.
            ## if label == "OPTION":
            ##     try:
            ##         _, name, value_str = fields
            ##     except ValueError:
            ##         raise FragmentFormatError("Expected 3 fields on OPTION line, not %d" % (len(fields),),
            ##                                   location)
            ##     if name not in options_dict:
            ##         print("options_dict", options_dict)
            ##         raise FragmentFormatError("Cannot set the new option %r" % (name,), location)
            ##     old_lineno, old_value_str = options_dict[name]
            ##     if old_value_str != value_str:
            ##         raise FragmentFormatError("Cannot modify option %r from %r to %r "
            ##                                   % (name, old_value_str, value_str), location)
            ##     continue

            raise FragmentFormatError("Unknown label %r" % (label,), location)

    finally:
        location.save(recno=get_recno(),
                      lineno=get_lineno(),
                      record=None)
        if close is not None:
            close()


class FileCache(object):
    def __init__(self, table, options):
        self.table = table
        self.options = options

    def get(self, name):
        return self.table.get(name)
            

def suggest_faster_json(reporter):
    loads = get_json_loads()
    if _json_module_name == "json":
        reporter.warning("Neither ujson nor cjson installed. Falling back to Python's slower built-in json decoder.")


def load_cache(filename, reporter):
    lc = filename.lower()
    if (lc.endswith(".fragments.gz") or lc.endswith(".fragments")):
        reporter = reporters.get_reporter(reporter)

        table = {}
        suggest_faster_json(reporter)
        with read_fragment_records(filename) as reader:
            for record in reporter.progress(reader, "Loading cache record"):
                table[record.id] = record

        return FileCache(table, reader.options)

    else:
        return open_fragdb(filename)

##### Write

def get_fragment_sort_key(frag):
    return (frag.num_cuts,
            frag.variable_symmetry_class, frag.variable_num_heavies, frag.variable_smiles,
            frag.constant_num_heavies, frag.constant_smiles, frag.constant_with_H_smiles,
            frag.attachment_order
            )


class FragmentWriter(object):
    def __init__(self, filename, outfile, options):
        self.filename = filename
        self._outfile = outfile
        self.options = options

    def close(self):
        self._outfile.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._outfile.close()

    def write_version(self):
        json.dump(("VERSION", "mmpdb-fragment/2"), self._outfile)
        self._outfile.write("\n")
        json.dump(("SOFTWARE", SOFTWARE), self._outfile)
        self._outfile.write("\n")

    def write_options(self, options):
        for k, v in sorted(options.to_text_settings()):
            if "\n" in k or "\r" in k or "\t" in k or " " in k:
                raise ValueError("Unsupported whitespace in key %r" % (k,))
            if "\n" in v or "\r" in v or "\t" in v:
                raise ValueError("Unsupported whitespace in %s value %r" % (k, v))
            json.dump(("OPTION", k, v), self._outfile)
            self._outfile.write("\n")

    def write_records(self, fragment_records):
        outfile = self._outfile
        
        for rec in fragment_records:
            if rec.errmsg:
                json.dump(("IGNORE", rec.id, rec.input_smiles, rec.errmsg), outfile)
                outfile.write("\n")
            else:
                fragment_fields = []
                record = ("RECORD", rec.id, rec.input_smiles, rec.num_normalized_heavies,
                          rec.normalized_smiles, fragment_fields)

                fragmentations = sorted(rec.fragments, key = get_fragment_sort_key)
                for frag in fragmentations:
                    fragment_fields.append((
                        frag.num_cuts, frag.enumeration_label,
                        frag.variable_num_heavies, frag.variable_symmetry_class, frag.variable_smiles,
                        frag.attachment_order, frag.constant_num_heavies, frag.constant_symmetry_class,
                        frag.constant_smiles, frag.constant_with_H_smiles,
                        ))
                json.dump(record, outfile)
                outfile.write("\n")

_wildcard_pat = re.compile(  re.escape("[*]")
                           + "|"
                           + re.escape("*"))


def relabel(smiles, order=None):
    input_smiles = smiles
    input_order = order

    if order is None:
        order = list(range(smiles.count("*")))
    else:
        assert not isinstance(order[0], int), ("Fix this for Python 3", order)
        order = [int(c) for c in order]

    def add_isotope_tag_to_wildcard(m):
        return "[*:%d]" % (order.pop(0)+1,)

    return _wildcard_pat.sub(add_isotope_tag_to_wildcard, smiles)

#### "fragdb" -- SQLite-based fragment file

SCHEMA_FILENAME = os.path.join(os.path.dirname(__file__), "fragment_schema.sql")
_schema_template = None    
def get_schema_template():
    global _schema_template
    if _schema_template is None:
        with open(SCHEMA_FILENAME) as infile:
            _schema_template = infile.read()
    return _schema_template

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
    db = sqlite3.connect(filename)
    c = db.cursor()
    # See if this is a database
    try:
        c.execute("SELECT COUNT(*) FROM record")
    except sqlite3.OperationalError as err:
        raise ValueError(f"{filename!r} does not appear to be a fragdb")
    try:
        c.execute("SELECT COUNT(*) FROM options")
    except sqlite3.OperationalError as err:
        raise ValueError(f"{filename!r} does not appear to be a fragdb")

    for (count,) in c:
        if count != 1:
            raise ValueError(f"{filename!r} does not appear to be a fragdb")
        break
    
    options = select_options(c)
    return FragDB(None, options, db, c)
    
#### options

_option_attrs = [field.name for field in dataclasses.fields(FragmentOptions)]
_option_cols = ["version"] + _option_attrs

_option_fields = ", ".join(_option_cols)
_option_qs = ",".join("?" * len(_option_cols))

_insert_options_sql = f"INSERT INTO options({_option_fields}) VALUES({_option_qs})"

def insert_options(c, options):
    values = [3] + [getattr(options, attr) for attr in _option_attrs]
    c.execute(_insert_options_sql, values)

def select_options(c):
    sql = f"SELECT {_option_fields} FROM options"
    for row in c.execute(sql):
        version, *values = row
        if version != 3:
            raise ValueError(f"Expected version 3 options, not version {version}")
        return FragmentOptions(*values)
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
        fragmentations = list(select_fragmentations_by_record_id(c, record_id)),
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
    record_c.execute(select_fragmentations_sql)
    for row in record_c:
        record_id, *record_values = row
        fragmentations = list(select_fragmentations_by_record_id(fragmentation_c, record_id))
        yield FragmentRecord(
            *record_values,
            fragmentations = fragmentations,
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
        self._record_id = 0 # Manage the record ids myself

    def close(self):
        db, c = self.db, self.c
        if db is not None:
            c.execute("CREATE INDEX fragmentation_on_record_id ON fragmentation(record_id)")
            c.execute("CREATE INDEX record_on_title ON record(title)")
            c.execute("CREATE INDEX error_record_on_title ON error_record(title)")
            self.db = self.c = None
            c.close()
            db.commit()  # don't rollback - keep partial writes too.
            db.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

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

def _get_fragdb_options(c):
    try:
        c.execute("SELECT count(*) FROM record")
    except sqlite3.OperationalError as err:
        raise ValueError(f"{filename!r} does not appear to be a fragdb file: {err}")

    fields = [
        "cut_smarts", "max_heavies", "max_rotatable_bonds",
        "method", "num_cuts", "rotatable_smarts",
        "salt_remover", "min_heavies_per_const_frag",
        ]
    field_str = ", ".join(f'"{name}"' for name in fields)
    try:
        c.execute(f"""
SELECT version, {field_str}
  FROM config
        """)
    except sqlite3.OperationalError as err:
        raise ValueError(f"{filename!r} does not appear to be a fragdb file: {err}")

    version, *field_values = next(c)
    if version != 3:
        raise ValueError(f"{filename!r} is a version {version} fragdb database. Only version 2 is supported.")

    kwargs = dict(zip(fields, field_values))
    return FragmentOptions(**kwargs)
                
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
    
### Dispatch to the correct writer

def open_fragment_writer(filename, options, format_hint=None):
    if format_hint is not None and format_hint not in ("fragdb", "fragments", "fragments.gz"):
        raise ValueError("Unsupported format_hint: %r" % (format_hint,))

    if format_hint is None:
        if filename is None:
            format_hint = "fragment"
        else:
            lc_filename = filename.lower()
            if ( lc_filename.endswith(".fragments.gz")
                or lc_filename.endswith(".fragments")):
                format_hint = "fragments"
            else:
                format_hint = "fragdb"
            
    if "fragments" in format_hint:
        outfile = fileio.open_output(filename, format_hint)
        writer = FragmentWriter(filename, outfile, options)
        writer.write_version()
        writer.write_options(options)
    else:
        if filename is None:
            filename = "input.fragdb"
        try:
            os.unlink(filename)
        except FileNotFoundError:
            pass
        db = sqlite3.connect(filename)
        c = db.cursor()
        writer = FragDBWriter(filename, db, c, options)
        init_fragdb(c, options)
        return writer

    # FragInfoWRiter and FragmentWriter but not FragDB
    return writer
