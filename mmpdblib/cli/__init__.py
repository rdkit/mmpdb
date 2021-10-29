from .. import __version__

import click

from .click_utils import ordered_group

from . import fragment
from . import smifrag
from . import index

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
        self.explain = reporter.explain
        
    
        
@ordered_group(epilog = epilog)
@click.option(
    "--quiet",
    "-q",
    is_flag=True,
    help="do not show progress or status information"
    )
@click.version_option(version = __version__)
@click.pass_context
def main(ctx, quiet):
    "Matched-molecular pair database loader"
    if ctx.obj is None:
        ctx.obj = CmdConfig(quiet)

main.add_command(fragment.fragment)
main.add_command(smifrag.smifrag)
main.add_command(index.index)

@main.command()
def spam():
    pass


## add_commands(cli_fragdb.fragdb_utils)

if __name__ == "__main__":
    main()
