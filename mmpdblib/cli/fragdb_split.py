
import os
import collections
import random

import click
from .click_utils import (
    command,
    die,
    positive_int,
    template_type,
    add_single_database_parameters,
    open_fragdb_from_options_or_exit,
    )


def get_constants_from_file(infile, has_header):
    constants = set()
    line_iter = iter(infile)

    if has_header:
        # Skip the first line, if it exists
        try:
            next(line_iter)
        except StopIteration:
            return constants
        
    for line in line_iter:
        if line[:1].isspace():
            continue
        fields = line.split(None, 1)
        assert fields, ("should not happen", line)
        smiles = fields[0]
        constants.add(smiles)
    return constants

def get_constant_counts(constants, fragdb, reporter):
    c = fragdb.cursor()
    constant_counts = []
    for constant in constants:
        c.execute(
            "SELECT count(*) FROM fragmentation WHERE constant_smiles = ?",
            (constant,))
        for (n,) in c:
            break
        else:
            raise AssertionError
        if n == 0:
            reporter.warning(f"Constant {constant!r} not in the database - skipping.")
        constant_counts.append((n, constant))
    return constant_counts
        
    
def get_constant_counts_from_db(db):
    c = db.cursor()
    c.execute("SELECT count(*), constant_smiles FROM fragmentation GROUP BY constant_smiles")
    return list(c)

def get_estimated_num_pairs(n):
    # +1 to include the possibility of matching to a hydrogen SMILES
    # and to keep from dumping everything into a single file.
    return n*(n-1)//2 + 1

def _largest_subset_sort_key(pair):
    tot_num_pairs, subset = pair
    return (-tot_num_pairs, len(subset), subset)

# first-fit to the subset with the smallest number of estimated pairs
def subset_by_max_files(constant_counts, num_files):
    assert num_files > 0, num_files
    import heapq
    heap = [(0, []) for i in range(num_files)]
    for count, constant in constant_counts:
        num_pairs = get_estimated_num_pairs(count)
        tot_num_pairs, subset = heapq.heappop(heap)
        subset.append(constant)
        heapq.heappush(heap, (tot_num_pairs + num_pairs, subset))

    heap.sort(key = _largest_subset_sort_key)
    return [subset for (_, subset) in heap]

# first-fit to the subset with the smallest number of estimated pairs
# but add a new subset if too full
def subset_by_max_pairs(constant_counts, max_pairs):
    assert max_pairs > 0, max_pairs
    import heapq
    heap = [(0, [])]
    
    for count, constant in constant_counts:
        num_pairs = get_estimated_num_pairs(count)

        tot_num_pairs, subset = heap[0]
        
        if (not tot_num_pairs) or (tot_num_pairs + num_pairs <= max_pairs):
            # Add if the smallest is empty, or if there's room.
            subset.append(constant)
            heapq.heapreplace(heap, (tot_num_pairs + num_pairs, subset))
        else:
            # Doesn't fit into the smallest available subset.
            # Need a new one.
            # (If num_pairs > max_pairs could append to a special 'full' list.)
            heapq.heappush(heap, (num_pairs, [constant]))
            
    heap.sort(key = _largest_subset_sort_key)
    return [subset for (_, subset) in heap]
    

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

# --num-files 10
# --max-estimated-pairs 100000
# 

SET_LIMIT = "mmpdblib.fragdb_split.limit"
def set_num_files(ctx, param, value):
    if value is not None:
        prev = ctx.meta.get(SET_LIMIT, None)
        if prev is None:
            ctx.meta[SET_LIMIT] = "num-files"
        elif prev != "num-files":
            raise click.UsageError("Cannot use both --num-files and --max-estimated-pairs")
    return value

def set_max_estimated_pairs(ctx, param, value):
    if value is not None:
        prev = ctx.meta.get(SET_LIMIT, None)
        if prev is None:
            ctx.meta[SET_LIMIT] = "max-estimated-pairs"
        elif prev != "max-estimated-pairs":
            raise click.UsageError("Cannot specify both --num-files and --set_max_estimated_pairs")
    return value
    

@command(
    epilog = fragdb_split_epilog,
    name = "fragdb_split",
    )

@click.option(
    "--num-files",
    "-n",
    type = positive_int(),
    default = None,
    callback = set_num_files,
    help = "maximum number of files to generate",
    )

@click.option(
    "--max-estimated-pairs",
    type = positive_int(),
    default = None,
    callback = set_max_estimated_pairs,
    help = "maximum number of estimated pairs per file",
    )


@click.option(
    "--constants",
    "-c",
    "constants_file",
    type = click.File("r"),
    help = "only export fragmentations containing the constants specified in the named file",
    )
@click.option(
    "--has-header / --no-header",
    "has_header",
    default = True,
    help = (
        "With --has-header (the default), skip the first line of the constants file. "
        "With --no-header, interpret the first line as a SMARTS"
        ),
    )

@click.option(
    "--template",
    "-t",
    default = "{prefix}.{i:04}.fragdb",
    type = template_type(),
    show_default = True,
    help = "template for the output filenames",
    )

@add_single_database_parameters()
@click.pass_obj
def fragdb_split(
        reporter,
        database_options,
        num_files,
        max_estimated_pairs,
        template,
        constants_file,
        has_header,
        ):
    """split a fragdb based on common constants

    DATABASE - a fragdb fragments database filename
    """
    import pathlib
    import shutil

    fragdb = open_fragdb_from_options_or_exit(database_options)

    # Get values used for template generation
    database_path = pathlib.Path(database_options.database)
    database = str(database_path)
    database_parent = database_path.parent
    database_stem = database_path.stem
    database_prefix = str(database_path.parent / database_path.stem)


    if constants_file is not None:
        with constants_file:
            constants = get_constants_from_file(constants_file, has_header)
        src = f"file {constants_file.name!r}"
        if not constants:
            reporter.report("No constants found in {src}. Exiting.")
            return
        constant_counts = get_constant_counts(constants, fragdb, reporter)
        if not constant_counts:
            reporter.report("No constants from {src} found in {database_options.database!r}. Exiting.")
            return
    else:
        constant_counts = get_constant_counts_from_db(fragdb)
        src = f"database {database_options.database!r}"
        if not constant_counts:
            reporter.report("No constants from {src}. Exiting.")
    
    reporter.report(f"Using {len(constant_counts)} constants from {src}.")

    # Sort from the most common to the least so we can use a first-fit algorithm.
    # (Decreasing count, increasing constant)
    constant_counts.sort(key = lambda pair: (-pair[0], pair[1]))
    
    # Assign to bins
    if num_files is None:
        if max_estimated_pairs is None:
            # Is this a good default?
            constant_subsets = subset_by_max_files(constant_counts, 10)
        else:
            constant_subsets = subset_by_max_pairs(constant_counts, max_estimated_pairs)
    else:
        constant_subsets = subset_by_max_files(constant_counts, num_files)
        
    
    # Save to files
    for i, subset in enumerate(constant_subsets):
        output_filename = template.format(
            parent = database_parent,
            stem = database_stem,
            prefix = database_prefix,
            sep = os.sep,
            i = i,
            )

        # Copy to the destination
        if not subset:
            continue
        reporter.report(f"Exporting {len(subset)} constants to {output_filename!r}")
        try:
            shutil.copy(database, output_filename)
        except IOError as err:
            die(f"Cannot copy {database!r} to {output_filename!r}: {err}")

        restrict_to_subset(output_filename, subset)
