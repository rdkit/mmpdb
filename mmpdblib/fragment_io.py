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

from . import __version__ as mmpdblib_version

from . import config
from . import fileio
from . import fragment_algorithm
from . import reporters
from .fragment_types import FragmentRecord, FragmentErrorRecord, FragmentFormatError
from ._compat import basestring

SOFTWARE = "mmpdb-" + mmpdblib_version

##### Read

parse_max_heavies_value = config.positive_int_or_none
parse_max_rotatable_bonds_value = config.positive_int_or_none

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
    }
    
def _get_options(line_reader, location):
    options_dict = {}
    options = config.FragmentOptions(**config.DEFAULT_FRAGMENT_OPTIONS.to_dict())
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
                    raise AssertionError # I don't know what caused this error

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
    reporter = reporters.get_reporter(reporter)

    table = {}
    suggest_faster_json(reporter)
    with read_fragment_records(filename) as reader:
        for record in reporter.progress(reader, "Loading cache record"):
            table[record.id] = record
        
    return FileCache(table, reader.options)
            
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
                        frag.attachment_order,
                        frag.constant_num_heavies, frag.constant_symmetry_class, frag.constant_smiles, frag.constant_with_H_smiles,
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
    
        
            
class FragInfoWriter(object):
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
        self._outfile.write("FORMAT mmpdb-fraginfo/2\n")
        self._outfile.write("SOFTWARE " + SOFTWARE + "\n")
        
        
    def write_options(self, options):
        for k, v in sorted(options.to_text_settings()):
            if "\n" in k or "\r" in k or "\t" in k or " " in k:
                raise ValueError("Unsupported whitespace in key %r" % (k,))
            if "\n" in v or "\r" in v or "\t" in v:
                raise ValueError("Unsupported whitespace in %s value %r" % (k, v))
            self._outfile.write("OPTION %s=%s\n" % (k, v))
            
    def write_records(self, fragment_records):
        outfile = self._outfile
        
        for rec in fragment_records:
            if rec.errmsg:
                outfile.write("IGNORE id=%r %r errmsg=%r\n"
                              % (rec.id, rec.input_smiles, rec.errmsg))
            else:
                outfile.write("RECORD id=%r %r #heavies=%d #fragmentations=%d\n"
                              % (rec.id, rec.input_smiles, rec.num_normalized_heavies, len(rec.fragments)))
                
                fragmentations = sorted(rec.fragments, key = get_fragment_sort_key)
                for frag in fragmentations:
                    reaction = "variable %s // constant %s" % (
                        relabel(frag.variable_smiles, frag.attachment_order),
                        relabel(frag.constant_smiles))
                                              
                    outfile.write(" FRAG #cuts=%d enum_label=%s %s\n"
                                  "   variable: #heavies=%d symm_class=%s %s attachment_order=%s\n"
                                  "   constant: #heavies=%d symm_class=%s %s H-smiles=%s\n"
                                  % (frag.num_cuts, frag.enumeration_label, reaction,
                                     frag.variable_num_heavies, frag.variable_symmetry_class,
                                     frag.variable_smiles, frag.attachment_order,
                                     frag.constant_num_heavies, frag.constant_symmetry_class,
                                     frag.constant_smiles, frag.constant_with_H_smiles
                                     ))

def open_fragment_writer(filename, options, format_hint=None):
    if format_hint is not None and format_hint not in ("fragments", "fragments.gz", "fraginfo", "fraginfo.gz"):
        raise ValueError("Unsupported format_hint: %r" % (format_hint,))
    outfile = fileio.open_output(filename, format_hint)

    if format_hint is None:
        if filename is None:
            format_hint = "fragment"
        else:
            lc_filename = filename.lower()
            if (   lc_filename.endswith(".fraginfo.gz")
                or lc_filename.endswith(".fraginfo")):
                format_hint = "fraginfo"
            else:
                format_hint = "fragment"
            
    if "fraginfo" in format_hint:
        writer = FragInfoWriter(filename, outfile, options)
    else:
        writer = FragmentWriter(filename, outfile, options)

    writer.write_version()
    writer.write_options(options)
    return writer
