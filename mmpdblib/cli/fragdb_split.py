
import os
import collections
import random

import click
from .click_utils import (
    command,
    die,
    positive_int_or_none,
    add_single_database_parameters,
    open_fragdb_from_options_or_exit,
    )


class template_type(click.ParamType):
    name = "STR"
    def convert(self, value, param, ctx):
        if value is None:
            return value

        reference_d = {
            "prefix": "blah",
            "parent": "blah",
            "stem": "blah",
            "sep": os.sep,
            "i": 1,
            }
        try:
            value.format(**reference_d)
        except ValueError as err:
            self.fail(f"Unable to evaluate: {err}", param, ctx)
        except KeyError as err:
            self.fail(f"Unsupported field: {err}", param, ctx)

        return value


def get_constants_from_file(infile):
    constants = set()
    for line in infile:
        if line[:1].isspace():
            continue
        fields = line.split(None, 1)
        assert fields, ("should not happen", line)
        smiles = fields[0]
        constants.add(smiles)
    return constants

def get_constants_from_db(db):
    c = db.cursor()
    c.execute("SELECT DISTINCT constant_smiles FROM fragmentation")
    constants = set(s for (s,) in c)
    return constants


def restrict_to_subset(output_filename, constant_smiles):
    import sqlite3

    db = sqlite3.connect(output_filename)
    c = db.cursor()
    
    # Create a new table containing the constant smiles
    c.execute("""
CREATE TEMPORARY TABLE required_constants (
  constant_smiles TEXT
)""")
    c.executemany("""
INSERT INTO required_constants (constant_smiles) VALUES (?)
""", ((smiles,) for smiles in constant_smiles)
     )
    ## # Index
    ## c.execute("CREATE INDEX required_constants_idx ON required_constants(constant_smiles)")

    # Remove fragmentations which don't have the required constant
    c.execute("""
DELETE FROM fragmentation
      WHERE constant_smiles NOT IN (
           SELECT constant_smiles FROM required_constants
         )
""")
    c.execute("COMMIT")
    

        
fragdb_split_epilog = """

The constants file is a line-oriented format with one constant
fragment SMILES per line. The SMILES must end with a whitepace or
a newline. Any characters after the whitespace are ignored.

The template is used to generate the output filenames, formatted
using [Format String Syntax](https://docs.python.org/3/library/string.html#formatstrings).
The available fields are:

\b
  {prefix} - the DATABASE path without the final extension
  {parent} - the DATABASE parent directory, or '.'
  {stem} - the DATABASE filename without the directory or final extension
  {sep} - the filesystem path seperator (eg, '/')
  {i} - an integer value 0 <= i < n

For example, if the DATABASE name is '/abc/xyz.fragdb', on a macOS
system, then:

\b
  {prefix} = '/abc/xyz'
  {parent} = '/abc'
  {stem} = 'xyz'
  {sep} = '/'


"""

@command(
    epilog = fragdb_split_epilog,
    name = "fragdb_split",
    )

@click.option(
    "-n",
    type = positive_int_or_none(),
    default = 10,
    help = "maximum number of files to generate (default: 10)",
    )


@click.option(
    "--constants",
    "constants_file",
    type = click.File("r"),
    help = "only export fragmentations containing the constants specified in the named file",
    )

@click.option(
    "--template",
    default = "{prefix}.{i:04}.fragdb",
    type = template_type(),
    show_default = True,
    help = "template for the output filenames",
    )

@click.option(
    "--seed",
    type = click.INT,
    help = "seed used to randomize how the constants are assigned to files",
    )

@add_single_database_parameters()
@click.pass_obj
def fragdb_split(
        reporter,
        database_options,
        n,
        template,
        seed,
        constants_file,
        ):
    """split a fragdb based on common constants

    DATABASE - a fragdb fragments database filename
    """
    import pathlib
    import shutil

    fragdb = open_fragdb_from_options_or_exit(database_options)

    database_path = pathlib.Path(database_options.database)
    database = str(database_path)
    database_parent = database_path.parent
    database_stem = database_path.stem
    database_prefix = str(database_path.parent / database_path.stem)
    
    if constants_file is not None:
        with constants_file:
            constants = get_constants_from_file(constants_file)
        src = f"file {constants_file.name!r}"
    else:
        constants = get_constants_from_db(fragdb)
        src = f"database {database_options.database!r}"
        
    if not constants:
        reporter.report("No constants found in {src}. Exiting with no output.")
        return
    
    reporter.report(f"Using {len(constants)} constants from {src}.")

    if n == "none":
        n = len(constants)

    # Randomize the contants
    constants = list(constants)
    random.Random(seed).shuffle(constants)

    # Assign to bins
    constant_bins = collections.defaultdict(set)
    for i, constant in enumerate(constants):
        constant_bins[i % n].add(constant)

    for i, subset in sorted(constant_bins.items()):
        output_filename = template.format(
            parent = database_parent,
            stem = database_stem,
            prefix = database_prefix,
            sep = os.sep,
            i = i,
            )

        # Copy to the destination
        reporter.report(f"Exporting {len(subset)} constants to {output_filename!r}")
        try:
            shutil.copy(database, output_filename)
        except IOError as err:
            die(f"Cannot copy {database!r} to {output_filename!r}: {err}")

        restrict_to_subset(output_filename, subset)

