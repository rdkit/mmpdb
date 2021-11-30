"Implement the 'rgroup2smarts' command"

# mmpdb - matched molecular pair database generation and analysis
#
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
#    * Neither the name of Andrew Dalke Scientific. nor the names of
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

# Convert a R-group SMILES into a SMARTS pattern.
# Must be rooted at "*", with one single bond.

# The algorithm is:
#  - figure out the valence and hydrogen count for each atom
#  - convert the molecule into an isomeric SMILES, with
#      hydrogens on the atom terms
#  - convert the SMILES into a SMARTS:
#     - use a regexp to find the atom terms
#     - insert with the valence

# For examples:
#   *-Cl -> *!@[Clv1]
#   *-O  -> *!@[OH1v2]

import sys

import click
from .click_utils import command, die

rgroup2smarts_epilog = """

This command is primarily meant for users to see how the `mmpdb
fragment` parameters `--cut-rgroup` and `--cut-rgroup-file` work.

A fragment file contains one fragment SMILES per line. Each fragment
SMILES must contain one and only one wildcard atom ("*"), which marks
the attachment point.

Blank lines and leading whitespace are not supported. The SMILES ends
at the first whitespace. Additional text on a line is ignored.

Each fragment is turned into a SMARTS pattern which matches that
fragment. By default the SMARTS patterns are converted into a
recursive SMARTS with all of the fragments. Use `--single` to output
the non-recursive SMARTS for each input SMILES.

Use `--check` to verify that the final SMARTS matches the input
fragments. Use `--cut-rgroup` to specify the SMILES fragments on the
command-line instead of from a file.
"""


@command(epilog=rgroup2smarts_epilog)
@click.option(
    "--cut-rgroup",
    metavar="SMILES",
    multiple=True,
    help="R-group SMILES to use",
)
@click.option(
    "--single",
    "-s",
    default=False,
    is_flag=True,
    help="Generate a SMARTS for each R-group SMILES (default: generate a single recursive SMARTS)",
)
@click.option(
    "--check",
    "-c",
    default=False,
    is_flag=True,
    help="Check that the SMARTS strings are valid (default: assume they are valid)",
)
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="Write conversion and check details to stderr",
)
@click.argument(
    "rgroup_filename",
    metavar="FILENAME",
    required=False,
)
@click.pass_obj
def rgroup2smarts(
    reporter,
    check,
    explain,
    cut_rgroup,
    rgroup_filename,
    single,
):
    """Convert an R-group file into a SMARTS which matches all of the SMILES

    FILENAME: file containing one or more R-group SMILES (use stdin if not specified)
    """
    from rdkit import Chem
    from .. import rgroup2smarts as _rgroup2smarts

    reporter.set_explain(explain)
    explain = reporter.explain

    close = None

    if cut_rgroup:
        if rgroup_filename is not None:
            die("Cannot specify both an R-group filename and a --cut-rgroup")
        location = _rgroup2smarts.ListLocation("--cut-rgroup SMILES")
        location.save(recno=1)
        explain("Using --cut-rgroup SMILES from the command-line")
        record_reader = _rgroup2smarts.iter_smiles_list(cut_rgroup, location)

    elif rgroup_filename is not None:
        explain(f"Reading R-group SMILES from {rgroup_filename!r}")
        location = _rgroup2smarts.FileLocation(rgroup_filename)
        try:
            f = open(rgroup_filename)
        except OSError as err:
            die(f"Cannot open input file: {err}")
        close = f.close
        record_reader = _rgroup2smarts.parse_rgroup_file(f, location)

    else:
        explain("Reading R-group SMILES from <stdin>")
        location = _rgroup2smarts.FileLocation("<stdin>")
        record_reader = _rgroup2smarts.parse_rgroup_file(sys.stdin, location)

    if check:
        all_mols = []
    else:
        all_mols = None

    outfile = sys.stdout

    iter_smarts = _rgroup2smarts.iter_smiles_as_smarts(record_reader, location, explain, all_mols)

    all_smarts = None

    try:
        if single:
            for smarts in iter_smarts:
                outfile.write(smarts + "\n")
        else:
            all_smarts = []
            for smarts in iter_smarts:
                assert smarts.startswith("*-!@"), (smarts, location)
                all_smarts.append(smarts)
            if not all_smarts:
                die(f"Cannot make a SMARTS: no SMILES strings found in {location.filename!r}")

    except _rgroup2smarts.ParseError as err:
        die(f"Cannot parse input file: {err}")
    except _rgroup2smarts.ConversionError as err:
        die(str(err))
    finally:
        if close is not None:
            close()

    if not single:
        smarts = _rgroup2smarts.make_recursive_smarts(all_smarts)

        try:
            if check:
                explain("Checking that the SMARTS matches all of the input molecules")
                all_pat = Chem.MolFromSmarts(smarts)
                if all_pat is None:
                    die(f"Cannot process final SMARTS: {smarts!r}")

                for i, (mol, where, smiles) in enumerate(all_mols):
                    if not mol.HasSubstructMatch(all_pat):
                        die(f"final SMARTS does not match SMILES from {where} ({smiles!r})")
                    explain(f"checked #{i}")
        finally:
            outfile.write(smarts + "\n")

    outfile.flush()
