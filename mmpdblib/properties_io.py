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


class Properties(object):
    def __init__(self, id_header, id_column, property_names, property_columns,
                 property_table):
        self.id_header = id_header
        self.id_column = id_column
        self.property_names = property_names
        self.property_columns = property_columns
        self.property_table = property_table
        
    def get_ids(self):
        return self.id_column
    
    def get_property_values(self, id):
        return self.property_table[id]
        
    def iter_properties(self):
        for name, column in zip(self.property_names, self.property_columns):
            yield name, column


# If a line contains a tab then it's tab-delimited.  Otherwise it's
# whitespace delimited.  (Originally it was whitespace delimited. Then
# I decided to support identifiers which contained a space in
# them. This seemed like a hacky-but-good-enough solution.)
def _split(line):
    if "\t" in line:
        return line.rstrip("\n").split("\t")
    return line.split()


def load_properties(properties_file, reporter):
    try:
        header_line = next(properties_file)
    except StopIteration:
        raise ValueError("First line of the properties file must contain the header")
    header_names = _split(header_line)
    if not header_names:
        raise ValueError("The properties file must contain at least one column name, for the id")
    if header_names[0] not in ("id", "ID", "Name", "name"):
        reporter.warning("the identifier column in the properties file (column 1) has "
                         "a header of %r; should be 'id', 'ID', 'Name', or 'name'" % (header_names[0],))

    seen = set()
    for header_name in header_names:
        if header_name in seen:
            raise ValueError(
                "Duplicate header %r found. A property name may not be listed more than once."
                % (header_name,))
        seen.add(header_name)
        
    n = len(header_names)

    id_column = []
    property_table = {}
    property_rows = []
    for lineno, line in enumerate(properties_file, 2):
        fields = _split(line)
        if len(fields) != n:
            raise ValueError("Line %d has %d fields but the header has %d"
                             % (lineno, len(fields), n))
        float_fields = []
        try:
            for field in fields[1:]:
                if field == "*":
                    float_fields.append(None)
                else:
                    float_fields.append(float(field))
        except ValueError:
            raise ValueError("Line %d value %r cannot be converted to a float"
                             % (lineno, field))

        id = fields[0]
        id_column.append(id)
        property_table[id] = float_fields
        property_rows.append(float_fields)

    if property_rows:
        property_columns = list(zip(*property_rows))
    else:
        property_columns = [[] for _ in header_names]
    return Properties(header_names[0], id_column, header_names[1:], property_columns,
                      property_table)
