# This package is a modified version of file readers in chemfp.

# Copyright (c) 2010-2016 Andrew Dalke Scientific, AB (Sweden)
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# The modifications to this file are covered under the mmpdb license.

from __future__ import print_function, absolute_import

import sys
import gzip
import io

from ._compat import basestring, open_universal, io_wrapper
from ._compat import binary_stdin, binary_stdout

class FileFormatError(ValueError):
    pass

class Outfile(object):
    def __init__(self, name, outfile, close):
        self.name = name
        self._outfile = outfile
        self._close = close

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if self._close is not None:
            close = self._close
            self._close = None
            close()

    def write(self, data):
        self._outfile.write(data)

    def writelines(self, lines):
        self._outfile.writelines(lines)
            

def open_input(filename):
    if filename is None:
        return sys.stdin
    
    if filename.endswith(".gz"):
        return io_wrapper(gzip.open(filename))
    return open_universal(filename)
        
def open_output(filename, format_hint):
    if format_hint is None:
        if isinstance(filename, basestring):
            is_compressed = filename.lower().endswith(".gz")
        else:
            is_compressed = False
    else:
        is_compressed = format_hint.lower().endswith(".gz")
        
    if filename is None:
        if is_compressed:
            outfile = gzip.GzipFile(fileobj=binary_stdout, mode="w")
            outfile = io_wrapper(outfile)
            return Outfile("<stdout>", outfile, outfile.close)
        else:
            return Outfile("<stdout>", sys.stdout, None)

    if is_compressed:
        outfile = io_wrapper(gzip.open(filename, "w"))
    else:
        outfile = open(filename, "w")
    return Outfile(filename, outfile, outfile.close)

## The Location object from chemfp

_where_template = {
    (0, 0, 0): "<unknown position>",
    (0, 0, 1): "record #%(recno)d",
    (0, 1, 0): "line %(lineno)d",
    (0, 1, 1): "line %(lineno)d, record #%(recno)d",
    (1, 0, 0): "file %(filename)r",
    (1, 0, 1): "file %(filename)r, record #%(recno)d",
    (1, 1, 0): "file %(filename)r, line %(lineno)d",
    (1, 1, 1): "file %(filename)r, line %(lineno)d, record #%(recno)d",
    }

class Location(object):
    """Get location and other internal reader and writer state information

    A Location instance gives a way to access information like
    the current record number, line number, and molecule object.::
    
      >>> with chemfp.read_molecule_fingerprints("RDKit-MACCS166",
      ...                        "ChEBI_lite.sdf.gz", id_tag="ChEBI ID") as reader:
      ...   for id, fp in reader:
      ...     if id == "CHEBI:3499":
      ...         print("Record starts at line", reader.location.lineno)
      ...         print("Record byte range:", reader.location.offsets)
      ...         print("Number of atoms:", reader.location.mol.GetNumAtoms())
      ...         break
      ... 
      [08:18:12]  S group MUL ignored on line 103
      Record starts at line 3599
      Record byte range: (138171, 141791)
      Number of atoms: 36

    The supported properties are:

      * filename - a string describing the source or destination
      * lineno - the line number for the start of the file
      * mol - the toolkit molecule for the current record
      * offsets - the (start, end) byte positions for the current record
      * output_recno - the number of records written successfully
      * recno - the current record number
      * record - the record as a text string
      * record_format - the record format, like "sdf" or "can"
       

    Most of the readers and writers do not support all of the properties.
    Unsupported properties return a None. The *filename* is a read/write
    attribute and the other attributes are read-only.
    
    If you don't pass a location to the readers and writers then they will
    create a new one based on the source or destination, respectively.
    You can also pass in your own Location, created as ``Location(filename)``
    if you have an actual filename, or ``Location.from_source(source)`` or
    ``Location.from_destination(destination)`` if you have a more generic
    source or destination.

    """
    _get_recno = None
    _get_output_recno = None
    _get_lineno = None
    _get_offsets = None
    _get_mol = None
    _get_record = None
    _record_format = None
    
    def __init__(self, filename=None):
        """Use *filename* as the location's filename"""
        self.filename = filename

    def __repr__(self):
        """Return a string like 'Location("<stdout>")'"""
        return "Location(%r)" % (self.filename,)

    def where(self):
        """Return a human readable description about the current reader or writer state.

        The description will contain the filename, line number, record
        number, and up to the first 40 characters of the first line of
        the record, if those properties are available.
        """
        filename = self.filename
        lineno = self.lineno
        recno = self.recno

        template = _where_template[ (filename is not None, lineno is not None, recno is not None) ]
        s = template % {"filename": filename, "lineno": lineno, "recno": recno}

        first_line = self.first_line

        if first_line:  # Don't show None and don't show blank lines
            if len(first_line) > 40:
                t = repr(first_line[:40])
                t = t[:-1] + " ..." + t[-1]
                s += ": first line starts %s" % (t,)
            else:
                s += ": first line is %r" % (first_line,)
        return s

    @classmethod
    def from_source(cls, source):
        """Create a Location instance based on the source

        If *source* is a string then it's used as the filename.
        If *source* is None then the location filename is "<stdin>".
        If *source* is a file object then its ``name`` attribute
        is used as the filename, or None if there is no attribute.
        """
        if source is None:
            return cls("<stdin>")
        if isinstance(source, basestring):
            return cls(source)
        return cls(getattr(source, "name", None))

    @classmethod
    def from_destination(cls, destination):
        """Create a Location instance based on the destination
        
        If *destination* is a string then it's used as the filename.
        If *destination* is None then the location filename is "<stdout>".
        If *destination* is a file object then its ``name`` attribute
        is used as the filename, or None if there is no attribute.
        """
        if destination is None:
            return cls("<stdout>")
        if isinstance(destination, basestring):
            return cls(destination)
        return cls(getattr(destination, "name", None))


    def clear_registry(self):
        """Part of the internal API, and subject to change."""
        self._get_recno = None
        self._get_output_recno = None
        self._get_lineno = None
        self._get_offsets = None
        self._get_mol = None
        self._get_record = None

    def register(self, **kwargs):
        """Part of the internal API, and subject to change."""
        for k, v in kwargs.items():
            if k in ("get_recno", "get_output_recno", "get_lineno", "get_offsets", "get_mol", "get_record"):
                setattr(self, "_" + k, v)
            else:
                raise KeyError(k)

    def save(self, **kwargs):
        """Part of the internal API, and subject to change."""
        for k, v in kwargs.items():
            if k in ("recno", "output_recno", "lineno", "offsets", "mol", "record"):
                def recall_value(value=v):
                    return value
                setattr(self, "_get_" + k, recall_value)
            elif k == "record_format":
                self._record_format = v
            else:
                raise KeyError(k)

    def get_registry(self):
        """Part of the internal API, and subject to change."""
        return {"get_recno": self._get_recno,
                "get_output_recno": self._get_output_recno,
                "get_lineno": self._get_lineno,
                "get_offsets": self._get_offsets,
                "get_mol": self._get_mol,
                "get_record": self._get_record,
                }

    @property
    def recno(self):
        """The current record number

        For writers this is the number of records sent to
        the writer, and output_recno is the number of records
        sucessfully written to the file or string.
        """
        _get_recno = self._get_recno
        if _get_recno is None:
            return None
        return _get_recno()

    @property
    def output_recno(self):
        """The number of records actually written to the file or string.

        The value ``recno - output_recno`` is the number of records
        sent to the writer but which had an error and could not be
        written to the output.
        """
        _get_output_recno = self._get_output_recno
        if _get_output_recno is None:
            return None
        return _get_output_recno()

    @property
    def lineno(self):
        """The current line number, starting from 1"""
        _get_lineno = self._get_lineno
        if _get_lineno is None:
            return None
        return _get_lineno()

    @property
    def offsets(self):
        """The (start, end) byte offsets, starting from 0

        *start* is the record start byte position and *end* is
        one byte past the last byte of the record.
        """
        _get_offsets = self._get_offsets
        if _get_offsets is None:
            return None
        return _get_offsets()

    @property
    def mol(self):
        """The molecule object for the current record"""
        _get_mol = self._get_mol
        if _get_mol is None:
            return None
        return _get_mol()

    @property
    def record(self):
        """The current record as an uncompressed text string"""
        _get_record = self._get_record
        if _get_record is None:
            return None
        return _get_record()

    @property
    def record_format(self):
        """The record format name"""
        return self._record_format

    @property
    def first_line(self):
        """The first line of the current record"""
        _get_record = self._get_record
        if _get_record is None:
            return None
        record = _get_record()
        if record is None:
            return None
        first_line, _, _ = record.partition("\n")
        return first_line.rstrip("\r")
        

#####
    
class SmilesReader(object):
    def __init__(self, reader, location):
        self.reader = reader
        self.location = location

    def close(self):
        self.reader = None

    def __iter__(self):
        if self.reader is None:
            raise ValueError("I/O operation on closed file")
        return self.reader

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.reader = None

        
def _read_smiles_file(infile, close, first_lineno, delimiter_message,
                      split_delimiter, location):
    lineno = 0
    line = None
    def get_recno():
        return lineno + 1
    def get_lineno():
        return lineno + first_lineno
    def get_record():
        return line

    location.register(
        get_recno = get_recno,
        get_lineno = get_lineno,
        get_record = get_record,
        )

    yield "ready"
    try:
        for lineno, line in enumerate(infile):
            terms = split_delimiter(line)
            if len(terms) < 2:
                raise FileFormatError(
                    "%s, %s" % (delimiter_message, location.where()))
            yield terms
    finally:
        location.save(
            recno = get_recno(),
            lineno = get_lineno(),
            record = None,
            )
        if close is not None:
            close()

def _split_whitespace(line):
    return line.split()

def _read_whitespace(infile, close, first_lineno, location):
    return _read_smiles_file(
        infile, close, first_lineno,
        "must contain at least two whitespace-delimited fields", _split_whitespace, location)


def _split_tab(line):
    if line[-1:] == "\n":
        line = line[:-1]
    return line.split("\t")

def _read_tab(infile, close, first_lineno, location):
    return _read_smiles_file(
        infile, close, first_lineno,
        "must contain at least two tab-delimited fields", _split_tab, location)

def _split_space(line):
    if line[-1:] == "\n":
        line = line[:-1]
    return line.split(" ")

def _read_space(infile, close, first_lineno, location):
    return _read_smiles_file(
        infile, close, first_lineno,
        "must contain at least two space-delimited fields", _split_space, location)

def _split_comma(line):
    if line[-1:] == "\n":
        line = line[:-1]
    return line.split(",")

def _read_comma(infile, close, first_lineno, location):
    return _read_smiles_file(
        infile, close, first_lineno,
        "must contain at least two comma-delimited fields", _split_comma, location)


def _split_to_eol(line):
    if line[-1:] == "\n":
        line = line[:-1]
    return line.split(None, 1)

def _read_to_eol(infile, close, first_lineno, location):
    return _read_smiles_file(
        infile, close, first_lineno,
        "must contain a whitespace to delimit the to-eol fields", _split_to_eol, location)


_delimiter_readers = {
    "whitespace": _read_whitespace,
    "tab": _read_tab,
    "space": _read_space,
    "to-eol": _read_to_eol,
    "comma": _read_comma,
    None: _read_whitespace,
    "native": _read_whitespace,
    }

def read_smiles_file(filename, format=None, delimiter="whitespace", has_header=False):
    if format is None:
        if filename is not None and filename.lower().endswith(".gz"):
            format = "smi.gz"
        else:
            format = "smi"
    elif format not in ("smi", "smi.gz"):
        if format is None:
            format = "smi"
        else:
            raise ValueError("Unsupported format: %r" % (format,))
    try:
        reader_function = _delimiter_readers[delimiter]
    except KeyError:
        raise ValueError("Unsupported delimiter: %r" % (delimiter,))

    if filename is None:
        if format == "smi.gz":
            infile = io_wrapper(gzip.GzipFile(fileobj=binary_stdin))
            close = None
        else:
            infile = sys.stdin
            close = None
        name = "<stdin>"
            
    else:
        if format == "smi.gz":
            infile = io_wrapper(gzip.GzipFile(filename))
            close = infile.close
        else:
            infile = open_universal(filename)
            close = infile.close
        name = filename

    location = Location.from_source(name)
    location.save(record_format="smi")
        
    first_lineno = 1
    if has_header:
        if infile.readline():
            first_lineno += 1

    reader = reader_function(infile, close, first_lineno, location)

    # initialize the reader by processing up to the first yield
    assert next(reader) == "ready"
    
    return SmilesReader(reader, location)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = "fileio test code")
    parser.add_argument("-i", "--in", dest="format")
    parser.add_argument("--has-header", action="store_true", default=False)
    parser.add_argument("--delimiter", default=None)
    parser.add_argument("filename", default=None, nargs="?")
    
    import sys
    args = parser.parse_args()
    with read_smiles_file(args.filename, args.format, args.delimiter, args.has_header) as reader:
        for x in reader:
            print(reader.location.filename, reader.location.lineno, reader.location.recno, x)
        
