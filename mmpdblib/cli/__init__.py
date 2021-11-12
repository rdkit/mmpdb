"implement the mmpdb command-line API"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2021, Andrew Dalke Scientific AB
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

from .. import __version__

import sys
import click

from .click_utils import ordered_group

from . import (
    fragment,
    smifrag,
    index,
    predict,
    transform,
    rgroup2smarts,
    list_,
    rulecat,
    loadprops,
    smicat,
    propcat,
    proprulecat,
    smi_split,
    fragdb_list,
    fragdb_constants,
    fragdb_partition,
    fragdb_merge,
    merge,
    drop_index,
    create_index,
    help_,
)


def add_commands(group):
    for name, cmd in group.commands.items():
        main.add_command(cmd, name=name)


epilog = """\
The 'mmpdb' program implements a set of subcommands to work with
a matched molecular pair database. The commands are roughly
categorized as "analysis" and "admin" commands.

The analysis commands fragment a SMILES database, index the fragments
into matched molecular pairs into a local SQLite database, import
molecular properties into a database, and searches the database for
possible transformations or predictions.

The admin commands are used to administer a database. They include
include ways to list the available data sets, dump the data as
a SMILES file or CSV file, update new properties, and
re-aggregate the rule dataset should any property values change.

For a short description of how to generate and use a dataset:

  % mmpdb help-analysis

For a short description of how to adminster a database:

  % mmpdb help-admin

The "help-*-format" commands (like "help-property-format") give
more details about the given format.

In addition, pass the "--help" option to a given command to see
the full list of options for the command.
"""


def explain(msg, *args):
    full_msg = (msg % args) + "\n"
    sys.stderr.write(full_msg)
    sys.stderr.flush()


def get_explain(use_explain, reporter=None):
    if use_explain:
        if reporter is None:
            return explain
        else:
            return reporter.explain

    from ..reporters import no_explain

    return no_explain


class CmdConfig:
    def __init__(self, quiet):
        self.quiet = quiet
        from .. import reporters

        if quiet:
            reporter = reporters.get_reporter("quiet")
        else:
            reporter = reporters.get_reporter("verbose")
        self.reporter = reporter

        self.report = reporter.report
        self.warning = reporter.warning
        self.progress = reporter.progress
        self.update = reporter.update
        self.explain = get_explain(False, self.reporter)

    def set_explain(self, use_explain):
        self.explain = get_explain(use_explain, self.reporter)


@ordered_group(epilog=epilog)
@click.option("--quiet", "-q", is_flag=True, help="do not show progress or status information")
@click.version_option(version=__version__)
@click.pass_context
def main(ctx, quiet):
    "Matched-molecular pair database loader"
    if ctx.obj is None:
        ctx.obj = CmdConfig(quiet)


main.add_command(fragment.fragment)
main.add_command(smifrag.smifrag)
main.add_command(index.index)
main.add_command(transform.transform)
main.add_command(predict.predict)
main.add_command(rgroup2smarts.rgroup2smarts)
main.add_command(list_.list_)
main.add_command(rulecat.rulecat)
main.add_command(loadprops.loadprops)
main.add_command(smicat.smicat)
main.add_command(propcat.propcat)
main.add_command(proprulecat.proprulecat)
main.add_command(smi_split.smi_split)
main.add_command(fragdb_list.fragdb_list)
main.add_command(fragdb_constants.fragdb_constants)
main.add_command(fragdb_partition.fragdb_partition)
main.add_command(fragdb_merge.fragdb_merge)
main.add_command(merge.merge)
main.add_command(drop_index.drop_index)
main.add_command(create_index.create_index)

add_commands(help_.help_group)
## add_commands(cli_fragdb.fragdb_utils)

if __name__ == "__main__":
    main()
