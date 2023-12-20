"Utilities for working with SMILES file input reading"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2019, Andrew Dalke Scientific, AB
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

from .click_utils import set_click_attrs


class SmiInputOptions:
    def __init__(self, format, delimiter, has_header):
        self.format = format
        self.delimiter = delimiter
        self.has_header = has_header
        
    def read_smiles_file(self, filename):
        from .. import fileio
        return fileio.read_smiles_file(
            filename,
            self.format,
            self.delimiter,
            self.has_header,
            )
        

def add_input_options(command):
    def add_option(*args, **kwargs):
        click.option(*args, **kwargs)(command)

    add_option(
        "--in",
        "-i",
        "in_format",
        type=click.Choice(["smi", "smi.gz"]),
        help=(
            "Input structuture format (one of 'smi', 'smi.gz'). "
            "If not specified, use the filename extension or default to 'smi'."
        ),
    )

    add_option(
        "--delimiter",
        default="whitespace",
        type=click.Choice(["whitespace", "to-eol", "comma", "tab", "space", "native"]),
        help=(
            "SMILES file delimiter style (one of 'whitespace' (default), 'to-eol', "
            "'comma', 'tab', 'space', or 'native')"
        ),
    )

    add_option(
        "--has-header",
        is_flag=True,
        default=False,
        help="Skip the first line, which is the header line",
    )

    # Wrap the command to convert the fragment option parameters
    # into a single object

    def make_input_options_wrapper(**kwargs):
        kwargs["input_options"] = SmiInputOptions(
            format=kwargs.pop("in_format"),
            delimiter=kwargs.pop("delimiter"),
            has_header=kwargs.pop("has_header"),
        )
        return command(**kwargs)

    set_click_attrs(make_input_options_wrapper, command)

    return make_input_options_wrapper
