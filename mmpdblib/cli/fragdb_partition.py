
import collections
import contextlib
import heapq
import os
import pathlib
import random
import sqlite3

import click
from .click_utils import (
    command,
    die,
    positive_int,
    template_type,
    add_multiple_databases_parameters,
    open_fragdb_from_options_or_exit,
    )

from .. import fragment_db
from .. import schema

class ValidityChecker:
    def __init__(self, first_database):
        self.first_database = first_database
        self.seen_titles = {}
        self.options = None

    def check(self, database, c):
        if self.options is None:
            self.options = fragment_db.select_options(c)
        else:
            options = fragment_db.select_options(c)
            self.check_for_bad_options(self.first_database, self.options, database, options)
        
        self.check_for_duplicate_titles(self.seen_titles, database, c)

    def check_for_bad_options(self, first_database, first_options, database, options):
        first_d = first_options.to_dict()
        d = options.to_dict()
        if d == first_d:
            return

        # Figure out which values are different
        lines = [
            f"Cannot merge. The fragment options in {database!r} differ from {first_database!r}."
            ]
        for k in d:
            if d[k] != first_d[k]:
                lines.append(f"  {k}: {d[k]!r} != {first_d[k]!r}")
        die(*lines)
    

    def check_for_duplicate_titles(self, seen_titles, database, c):
        c.execute("SELECT title FROM record")
        for title, in c:
            if title in seen_titles:
                prev_field, prev_database = seen_titles[title]
                die(
                    f"Cannot merge. Duplicate record {title!r} found as record in {database!r} and "
                    f"{prev_field} in {prev_database!r}"
                    )
            seen_titles[title] = ("record", database)

        c.execute("SELECT title FROM error_record")
        for title, in c:
            if title in seen_titles:
                prev_field, prev_database = seen_titles[title]
                die(
                    f"Cannot merge. Duplicate record {title!r} found as error record in {database!r} and "
                    f"{prev_field} in {prev_database!r}"
                    )
            seen_titles[title] = ("error record", database)
        
        

def get_all_constant_counts(databases, reporter):
    num_databases = len(databases)
    assert num_databases > 0
    constant_counts = collections.defaultdict(int)

    checker = ValidityChecker(databases[0])
    
    for database in databases:
        with contextlib.closing(open_fragdb_from_options_or_exit(database)) as db:
            with contextlib.closing(db.cursor()) as c:
                checker.check(database, c)
                
                num_constant_smiles = 0
                num_fragmentations = 0
                c.execute(
                    "SELECT constant_smiles, COUNT(*) FROM fragmentation GROUP BY constant_smiles"
                    )
                for constant_smiles, n in c:
                    constant_counts[constant_smiles] += n
                    num_constant_smiles += 1
                    num_fragmentations += n
                    
                reporter.report(
                    f"Analyzed {database!r}: "
                    f"#constants: {num_constant_smiles} "
                    f"#fragmentations: {num_fragmentations}"
                    )
    if num_databases > 1:
        reporter.report(
            f"Analyzed {num_databases} databases. Found "
            f"#constants: {len(constant_counts)} "
            f"#fragmentations: {sum(constant_counts.values())}"
            )
    
    return list(constant_counts.items())
    
    
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



def get_specified_constant_counts(constants, databases, reporter):
    constant_counts = dict((constant, 0) for constant in constants)

    # Use an in-memory database to help with the selection.
    constants_db = sqlite3.connect(":memory:")
    constants_c = constants_db.cursor()
    constants_c.execute("""
CREATE TABLE constant (
  constant_smiles TEXT
)""")
    constants_c.executemany(
        "INSERT INTO constant (constant_smiles) VALUES (?)",
        [(constant,) for constant in constants]
        )
    constants_c.execute("CREATE UNIQUE INDEX constant_on_smiles ON constant(constant_smiles)")
    constants_c.execute("COMMIT")

    checker = ValidityChecker(databases[0])
    try:
        for database in databases:
            # Verify it's valid
            with open_fragdb_from_options_or_exit(database) as db:
                checker.check(database, db.cursor())

            constants_c.execute("ATTACH DATABASE ? AS fragdb", (database,))
            try:
                constants_c.execute("""
  SELECT fragmentation.constant_smiles, COUNT(*)
    FROM fragdb.fragmentation as fragmentation,
         constant
   WHERE constant.constant_smiles = fragmentation.constant_smiles
GROUP BY fragmentation.constant_smiles
""")
                num_constant_smiles = 0
                num_fragmentations = 0
                for constant_smiles, n in constants_c:
                    constant_counts[constant_smiles] += n
                    num_constant_smiles += 1
                    num_fragmentations += n
                
                reporter.report(
                    f"Counts from {database!r}: "
                    f"#constants: {num_constant_smiles} "
                    f"#fragmentations: {num_fragmentations}"
                    )
                
            finally:
                constants_c.execute("DETACH DATABASE fragdb")
                
    finally:
        constants_c.close()
        constants_db.close()


    num_databases = len(databases)
    if num_databases > 1:
        reporter.report(
            f"Counts from {num_databases} databases. Found "
            f"#constants: {len(constant_counts)} "
            f"#fragmentations: {sum(constant_counts.values())}"
            )

    final_counts = []
    for constant, count in constant_counts.items():
        if count == 0:
            reporter.warning(f"Constant {constant!r} not in the database - skipping.")
        else:
            final_counts.append((constant, count))
    return final_counts
        

def get_weight(n):
    # +1 to include the possibility of matching to a hydrogen SMILES
    # and to keep from dumping everything into a single file.
    return n*(n-1)//2 + 1

def _largest_subset_sort_key(pair):
    tot_weight, subset = pair
    return (-tot_weight, len(subset), subset)

# first-fit to the subset with the smallest weight
def subset_by_max_files(constant_counts, num_files):
    assert num_files > 0, num_files
    heap = [(0, []) for i in range(num_files)]
    for constant, count in constant_counts:
        weight = get_weight(count)
        tot_weight, subset = heapq.heappop(heap)
        subset.append(constant)
        heapq.heappush(heap, (tot_weight + weight, subset))

    heap.sort(key = _largest_subset_sort_key)
    return heap

# first-fit to the subset with the smallest weight
# but add a new subset if too full
def subset_by_max_weight(constant_counts, max_weight):
    assert max_weight > 0, max_weight
    heap = [(0, [])]
    
    for constant, count in constant_counts:
        weight = get_weight(count)

        tot_weight, subset = heap[0]
        
        if (not tot_weight) or (tot_weight + weight <= max_weight):
            # Add if the smallest is empty, or if there's room.
            subset.append(constant)
            heapq.heapreplace(heap, (tot_weight + weight, subset))
        else:
            # Doesn't fit into the smallest available subset.
            # Need a new one.
            # (If weight > max_weight could append to a special 'full' list.)
            heapq.heappush(heap, (weight, [constant]))
            
    heap.sort(key = _largest_subset_sort_key)
    return heap
    

def init_output_database(output_c, subset):
    # Create the schema
    schema._execute_sql(output_c, fragment_db.get_schema_template())
    
    # Create some in-memory/temporary tables
    schema._execute_sql(output_c, """
ATTACH DATABASE ":memory:" AS merge;

CREATE TABLE merge.required_constants (
  constant_smiles TEXT
);
    """)

    # Add the constants for this database subset
    output_c.executemany("""
INSERT INTO merge.required_constants (constant_smiles) VALUES (?)
""", ((smiles,) for smiles in subset)
     )
    
    # Index
    output_c.execute("CREATE INDEX merge.required_constants_idx ON required_constants(constant_smiles)")


def check_for_duplicates(output_c, input_filename, output_filename, reporter):
    # Check for duplicate record titles 
    output_c.execute("""
SELECT old_record.title
  FROM old.record AS old_record,
       record AS new_record
 WHERE old_record.title = new_record.title
 LIMIT 1
""")
    for (title,) in output_c:
        reporter.update("")
        die(
            f"ERROR! Cannot merge {input_filename!r}. "
            f"Record {title!r} already in {output_filename!r}."
            )
        
    # Check for duplicate error record titles
    output_c.execute("""
SELECT old_error_record.title
  FROM old.error_record AS old_error_record,
       error_record AS new_error_record
 WHERE old_error_record.title = new_error_record.title
 LIMIT 1
""")
    for (title,) in output_c:
        die(
            f"Cannot merge {input_filename!r} ",
            f"Error record {title!r} already in {output_filename!r}."
            )
    
    
def export_options(output_c):
    # Copy the options
    output_c.execute("INSERT INTO options SELECT * FROM old.options")
    
def export_subset(output_c):
    # Copy the erorr_record rows
    output_c.execute("""
INSERT INTO error_record (title, input_smiles, errmsg)
SELECT title, input_smiles, errmsg
  FROM old.error_record
""")
    
    # Copy the record rows
    output_c.execute("""
INSERT INTO record (title, input_smiles, num_normalized_heavies, normalized_smiles)
SELECT title, input_smiles, num_normalized_heavies, normalized_smiles
  FROM old.record
""")

    # Copy the relevant fragmentations
    output_c.execute("""
INSERT INTO fragmentation (
    record_id,
    num_cuts,
    enumeration_label,
    variable_num_heavies,
    variable_symmetry_class,
    variable_smiles,
    attachment_order,
    constant_num_heavies,
    constant_symmetry_class,
    constant_smiles,
    constant_with_H_smiles) 
  SELECT record_id,
         num_cuts,
         enumeration_label,
         variable_num_heavies,
         variable_symmetry_class,
         variable_smiles,
         attachment_order,
         constant_num_heavies,
         constant_symmetry_class,
         old_fragmentation.constant_smiles,
         constant_with_H_smiles
    FROM old.fragmentation AS old_fragmentation,
         merge.required_constants AS required_constants
   WHERE old_fragmentation.constant_smiles = required_constants.constant_smiles
   ;
""")
    

        
fragdb_partition_epilog = """

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
# --max-weight 100000
# 

SET_LIMIT = "mmpdblib.fragdb_partition.limit"
def set_num_files(ctx, param, value):
    if value is not None:
        prev = ctx.meta.get(SET_LIMIT, None)
        if prev is None:
            ctx.meta[SET_LIMIT] = "num-files"
        elif prev != "num-files":
            raise click.UsageError("Cannot use both --num-files and --max-weight")
    return value

def set_max_weight(ctx, param, value):
    if value is not None:
        prev = ctx.meta.get(SET_LIMIT, None)
        if prev is None:
            ctx.meta[SET_LIMIT] = "max-weight"
        elif prev != "max-weight":
            raise click.UsageError("Cannot specify both --num-files and --max-weight")
    return value
    

@command(
    epilog = fragdb_partition_epilog,
    name = "fragdb_partition",
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
    "--max-weight",
    type = positive_int(),
    default = None,
    callback = set_max_weight,
    help = (
        "maximum weight per file (weight = N*(N-1)/2+1 where N "
        "is the number of occurance of the constant"
        ),
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

@click.option(
    "--reference",
    metavar = "STR",
    default = None,
    help = "reference filename to use for the template",
    )

@add_multiple_databases_parameters()
@click.pass_obj
def fragdb_partition(
        reporter,
        databases_options,
        num_files,
        max_weight,
        template,
        reference,
        constants_file,
        has_header,
        ):
    """partition a fragdb based on common constants

    DATABASE - a fragdb fragments database filename
    """

    databases = databases_options.databases
    if not databases:
        raise click.BadArgumentUsage("must specify at least one fragment database")
    
    num_databases = len(databases)

    if reference is None:
        if num_databases == 1:
            reference = databases[0]
        else:
            reference = "partition.fragdb"

    # Get values used for template generation
    reference_path = pathlib.Path(reference)
    reference_str = str(reference_path)
    reference_parent = reference_path.parent
    reference_stem = reference_path.stem
    reference_prefix = str(reference_path.parent / reference_path.stem)    

    # Load the constant counts
    if constants_file is None:
        constant_counts = get_all_constant_counts(databases, reporter)
    else:
        with constants_file:
            constants = get_constants_from_file(constants_file, has_header)
        if not constants:
            reporter.report("No constants found in {constants_file.name!r}. Exiting.")
            return            
        
        constant_counts = get_specified_constant_counts(constants, databases, reporter)

    # Sort from the most common to the least so we can use a first-fit algorithm.
    # (Decreasing count, increasing constant)
    constant_counts.sort(key = lambda pair: (-pair[1], pair[0]))
    
    # Assign to bins
    if num_files is None:
        if max_weight is None:
            # Is 10 a good default?
            constant_subsets = subset_by_max_files(constant_counts, 10)
        else:
            constant_subsets = subset_by_max_weight(constant_counts, max_weight)
    else:
        constant_subsets = subset_by_max_files(constant_counts, num_files)
        
    
    # Save to files
    for subset_i, (total_weight, subset) in enumerate(constant_subsets):
        output_filename = template.format(
            parent = reference_parent,
            stem = reference_stem,
            prefix = reference_prefix,
            sep = os.sep,
            i = subset_i,
            )

        # Copy to the destination
        if not subset:
            continue
        reporter.report(
            f"Exporting {len(subset)} constants to {output_filename!r} "
            f"(#{subset_i+1}/{len(constant_subsets)}, weight: {total_weight})",
            )

        try:
            os.unlink(output_filename)
        except IOError:
            pass

        try:
            output_db = sqlite3.connect(output_filename)
        except sqlite3.OperationalError as err:
            die("Cannot create output file {output_filename!r}: {err}")

        output_c = output_db.cursor()
        
        init_output_database(output_c, subset)
        output_c.execute("COMMIT")
        
        for database_i, database in enumerate(databases, 1):

            try:
                output_c.execute("ATTACH DATABASE ? AS old", (database,))
            except sqlite3.OperationalError as err:
                die("Cannot attach file {database!t} to {output_filename!r}: {err}")

            try:
                if num_databases > 1:
                    reporter.update(f"Exporting constants from {database!r} (#{database_i}/{num_databases})")
                    
                #check_for_duplicates(output_c, database, output_filename, reporter)
                
                output_c.execute("BEGIN TRANSACTION")
                if database_i == 1:
                    export_options(output_c)
                export_subset(output_c)
            finally:
                reporter.update("")
                try:
                    output_c.execute("COMMIT")
                except sqlite3.OperationalError:
                    pass
                output_c.execute("DETACH DATABASE old")

        # And index (need the title indices for duplicate checks.)
        schema._execute_sql(output_c, fragment_db.get_fragment_create_index_sql())
                
        output_db.close()
