"Implement the 'propcat' command"

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

import click

from .click_utils import (
    command,
    die,
    add_single_database_parameters,
    add_multiple_properties,
    open_database_from_options_or_exit,
    get_property_names_or_error,
)

propcat_epilog = """

Write information about the properties for the compounds in DATABASE,
formatted as a property file. Use `mmpdb help-property-file` for
details about the property file format.

The output from this command is a tab-delimited CSV file where the
first column has the head "ID" and contains the compound identifier.
The other columns contain property information for each compound. The
column title is the property name.

By default there is one column for each property in the databases, and
the one row for each compound with at least one property. Use
'--property' to limit the output to a specific property, or use it
multiple times to specify multiple property names to output. Use
'--all' to list all of the compounds, even if the compound has none of
the specified properties.

The character "*" will be use if a listed compound is missing a given
property.

Examples:

1) Write all of the properties to stdout:

\b
  % mmpdb propcat CHEMBL_thrombin_Ki_IC50.mmpdb

2) Write the "MP" property to "MP.properties":

\b
  % mmpdb propcat csd.mmpdb --property MP -o MP.properties

3) Write the compound identifiers only to stdout:

\b
  % mmpdb propcat csd.mmpdb --no-properties --all

"""


@command(epilog=propcat_epilog)
@add_single_database_parameters()
@add_multiple_properties
@click.option(
    "--all",
    "show_all",
    is_flag=True,
    default=False,
    help="Include compounds which have no properties",
)
@click.option(
    "--output",
    "-o",
    "output_filename",
    metavar="FILENAME",
    help="Output filename (default is stdout)",
)
@click.pass_obj
def propcat(
    reporter,
    database_options,
    show_all,
    property_names,
    no_properties,
    output_filename,
):
    """Write the database properties to a properties file

    DATABASE: an mmpdb file
    """
    from .. import fileio

    db = open_database_from_options_or_exit(database_options)
    c = db.get_cursor()
    dataset = db.get_dataset()

    property_names = get_property_names_or_error(
        dataset,
        property_names=property_names,
        no_properties=no_properties,
    )

    property_values_list = []
    for property_name in property_names:
        property_name_id = dataset.get_property_name_id(property_name)
        property_values_list.append(dataset.get_property_values(property_name_id))

    try:
        outfile = fileio.open_output(output_filename, output_filename)
    except IOError as err:
        die(f"Cannot open --output: {err}")

    with outfile:
        print("ID", *property_names, sep="\t", file=outfile)
        for compound_row in dataset.iter_compounds(cursor=c):
            columns = [compound_row.public_id]
            is_empty = True
            for property_values in property_values_list:
                value = property_values.get(compound_row.id, None)
                if value is None:
                    columns.append("*")
                else:
                    columns.append(value)
                    is_empty = False
            if show_all or not is_empty:
                print(*columns, sep="\t", file=outfile)
