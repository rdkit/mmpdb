"Implement the 'smifrag' command"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2019-2021, Andrew Dalke Scientific, AB
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
    name_to_command_line,
)

from .. import smarts_aliases

from . import fragment_click
from .. import fragment_records
from .. import fragment_types


########
smifrag_epilog = """
\b
""" + smarts_aliases.get_epilog(
    "--cut-smarts", smarts_aliases.cut_smarts_aliases
)


@command(epilog=smifrag_epilog)
@fragment_click.add_fragment_options
@click.argument(
    "smiles",
    # help="SMILES string to fragment"
)
@click.pass_context
def smifrag(ctx, fragment_options, smiles):
    """Fragment a single SMILES string

    SMILES: the SMILES string of the structure to fragment

    Fragment a SMILES and print details about each variable and constant
    fragment and how they are connected.

    """

    reporter = ctx.obj

    try:
        fragment_filter = fragment_options.get_fragment_filter()
    except fragment_types.FragmentValueError as err:
        die(f"Error in command-line option {name_to_command_line(err.name)!r} " f"({err.value!r}): err.reason")

    record = fragment_records.make_fragment_record_from_smiles(
        smiles,
        fragment_filter,
        reporter=reporter,
    )
    if record.errmsg:
        die(f"Cannot process SMILES: {record.errmsg}")

    columns = [
        ["#cuts"],
        ["enum.label"],
        ["#heavies"],
        ["symm.class"],
        ["smiles"],
        ["order"],
        ["#heavies"],
        ["symm.class"],
        ["smiles"],
        ["with-H"],
    ]
    styles = [
        "center",
        "center",
        "right",
        "center",
        "left",
        "center",
        "right",
        "center",
        "left",
        "left",
    ]

    for frag in record.fragmentations:
        items = [
            str(frag.num_cuts),
            frag.enumeration_label,
            str(frag.variable_num_heavies),
            frag.variable_symmetry_class,
            frag.variable_smiles,
            frag.attachment_order,
            str(frag.constant_num_heavies),
            frag.constant_symmetry_class,
            frag.constant_smiles,
            frag.constant_with_H_smiles or "-",
        ]

        for (item, column) in zip(items, columns):
            column.append(str(item))

    sizes = []
    for style, column in zip(styles, columns):
        column_width = max(map(len, column))
        column[0] = column[0].ljust(column_width)
        sizes.append(column_width)
        if len(column) == 1:
            continue

        data_width = max(map(len, column[1:]))

        if style == "center":
            column[1:] = [s.rjust(data_width).center(column_width) for s in column[1:]]
        elif style == "left":
            column[1:] = [s.ljust(data_width) for s in column[1:]]
        ## elif style == "right-10":
        ##     spacer = " "*10
        ##     column[1:] = [spacer + s.rjust(data_width-10).center(column_width-10) for s in column[1:]]
        elif style == "right":
            column[1:] = [s.rjust(data_width).center(column_width) for s in column[1:]]
        else:
            raise AssertionError(style)

    first_line = (
        " " * sizes[0]
        + "   "
        + " " * sizes[1]
        + " |-"
        + "  variable  ".center(sizes[2] + sizes[3] + sizes[4] + 6, "-")
        + "-| "
        + " " * sizes[5]
        + " |-"
        + "  constant  ".center(sizes[6] + sizes[7] + sizes[8] + sizes[9] + 9, "-")
    )

    print(first_line)
    for lineno, fields in enumerate(zip(*columns)):
        print(*fields, sep=" | ")
        if lineno == 0:
            print(*["-" * len(s) for s in fields], sep="-+-")
