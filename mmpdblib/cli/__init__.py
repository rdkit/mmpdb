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

from .click_utils import FormatEpilog

# Map from command to the mmpdb.cli.{module_name}.{function_name} handler.
# This is used to load the appropriate command dynamically.

command_groups = [
    (
        "Matched molecular pair generation commands (see 'help-analysis')", [
            ("fragment", "fragment.fragment"),
            ("smifrag", "smifrag.smifrag"),
            ("index", "index.index"),
            ("predict", "predict.predict"),
            ("transform", "transform.transform"),
            ("rgroup2smarts", "rgroup2smarts.rgroup2smarts"),
            ("generate", "generate.generate"),
            ],
    ), (
        "Distributed generation commands (see 'help-distributed')", [
            ("smi_split", "smi_split.smi_split"),
            ("fragdb_constants", "fragdb_constants.fragdb_constants"),
            ("fragdb_partition", "fragdb_partition.fragdb_partition"),
            ("fragdb_merge", "fragdb_merge.fragdb_merge"),
            ("merge", "merge.merge"),
            ],
    ), (
        "Administration commands (see 'help-admin')", [
            ("fragdb_list", "fragdb_list.fragdb_list"),
            ("list", "list_.list_"),
            ("loadprops", "loadprops.loadprops"),
            ("smicat", "smicat.smicat"),
            ("rulecat", "rulecat.rulecat"),
            ("ruleenvcat", "ruleenvcat.ruleenvcat"),
            ("propcat", "propcat.propcat"),
            ("proprulecat", "proprulecat.proprulecat"),
            ("drop_index", "drop_index.drop_index"),
            ("create_index", "create_index.create_index"),
            ],
    ), (
        "Help commands", [
            ("help", "help_.help_"),
            ("help-analysis", "help_.help_analysis"),
            ("help-admin", "help_.help_admin"),
            ("help-distributed", "help_.help_distributed"),
            ("help-postgres", "help_.help_postgres"),
            ("help-smiles-format", "help_.help_smiles_format"),
            ("help-property-format", "help_.help_property_format"),
            ]
    )
    ]
_commands = {}
for (title, command_pairs) in command_groups:
    _commands.update(command_pair for command_pair in command_pairs)
del title, command_pairs


epilog = """\
The 'mmpdb' program implements a set of subcommands to work with
a matched molecular pair database. The commands are roughly
categorized as "analysis", "distributed", and "admin" commands.

The analysis commands fragment a SMILES database, index the fragments
into matched molecular pairs into a local SQLite database, import
molecular properties into a database, and searches the database for
possible transformations or predictions.

The distributed commands are used to parallelize matched molecular
pair generation in a distributed cluster environment. They include
ways to split a SMILES file into parts, for parallel fragmentation, to
merge the fragmentations together for later re-use as a cache, to
partition the fragmentations by constant for parallel indexing, and to
merge the generated indices for a final database.

The admin commands are used to administer a database. They include
include ways to list the available data sets, dump the data as
a SMILES file or CSV file, update new properties, and
re-aggregate the rule dataset should any property values change.

For a short description of how to generate and use a dataset:

  % mmpdb help-analysis

See "help-postgres" for examples using Postgres.

For a short description of how to parallize MMP generation:

  % mmpdb help-distributed

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
        

# The 'main' command uses the MultiCommand click group.
# This lets me use my own methods to resolve the available commands.
class MultiCommand(FormatEpilog, click.MultiCommand):
    def list_commands(self, ctx):
        return list(_commands)

    def get_command(self, ctx, name):
        # secret alias
        if name == "transmogrify":
            name = "generate"
        path = _commands.get(name, None)
        if path is None:
            # Return None to indicate the command is not available.
            return None
        # Import the submodule then get and return the correct function.
        import importlib
        module_name, func_name = path.split(".")
        mod = importlib.import_module("." + module_name, __name__)
        return getattr(mod, func_name)
    
    def format_commands(self, ctx, formatter):
        # Use my own command formatter so I can group them by theme.
        limit = formatter.width - 6 - max(len(cmd[0]) for cmd in _commands)

        old_value = formatter.indent_increment
        formatter.indent_increment = 1
        try:
            with formatter.section("Commands"):
                for title, command_pairs in command_groups:
                    rows = []
                    for subcommand, _ in command_pairs:
                        cmd = self.get_command(ctx, subcommand)
                        help = cmd.get_short_help_str(limit)
                        rows.append((subcommand, help))

                    assert rows
                    with formatter.section(title):
                        formatter.write_dl(rows)
        finally:
            formatter.indent_increment = old_value
        
@click.group(epilog=epilog, cls=MultiCommand)
@click.option("--quiet", "-q", is_flag=True, help="Do not show progress or status information")
@click.version_option(version=__version__)
@click.pass_context
def main(ctx, quiet):
    "Matched-molecular pair database loader"
    if ctx.obj is None:
        ctx.obj = CmdConfig(quiet)

if __name__ == "__main__":
    main()
