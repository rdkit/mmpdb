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
    nonnegative_int,
    GzipFile,
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
            f"Cannot partition. The fragment options in {database!r} differ from {first_database!r}."
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
        
def validate_databases(databases):
    checker = ValidityChecker(databases[0])
    for database in databases:
        with contextlib.closing(open_fragdb_from_options_or_exit(database)) as db:
            with contextlib.closing(db.cursor()) as c:
                checker.check(database, c)

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
    
    return constant_counts
    
    
def get_constant_counts_from_file(infile, has_header, include_count):
    constant_counts = {}
    line_iter = iter(infile)
    lineno = 1
    seen = {}

    filename = infile.name
    def oopsie(msg):
        die(f"Unable to process constants file: {msg}, {filename!r} line {lineno}")
    
    if has_header:
        # Skip the first line, if it exists
        try:
            next(line_iter)
        except StopIteration:
            return constants
        lineno += 1
        
    for lineno, line in enumerate(line_iter, lineno):
        fields = line.split()
        num_fields = len(fields)
        if include_count:
            if num_fields == 0:
                oopsie("Missing SMILES and count columns")
            if num_fields == 1:
                oopsie("Missing count column")
            smiles, count_str = fields[:2]
            try:
                count = int(count_str)
                if count < 1:
                    raise ValueError
            except ValueError as err:
                oopsie(f"Count must be a positive integer (not {count_str!r})")
                
        else:
            if num_fields == 0:
                oopsie("Missing SMILES column")
            smiles = fields[0]
            count = 1

        if smiles in seen:
            prev_lineno = seen[smiles]
            oopsie(f"SMILES is a duplicate from {prev_lineno}")
        seen[smiles] = lineno
        
        constant_counts[smiles] = count
    return constant_counts


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
    schema._execute_sql(output_c, """
-- Copy the error_record rows
INSERT INTO error_record (title, input_smiles, errmsg)
SELECT title, input_smiles, errmsg
  FROM old.error_record;
    
-- Copy the record rows
INSERT INTO record (title, input_smiles, num_normalized_heavies, normalized_smiles)
SELECT title, input_smiles, num_normalized_heavies, normalized_smiles
  FROM old.record;

-- Copy the relevant fragmentation with this constant
-- Need to map to the correct record id

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
  SELECT new_record.id,
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
         merge.required_constants AS required_constants,
         record AS new_record,
         old.record AS old_record
   WHERE old_fragmentation.constant_smiles = required_constants.constant_smiles
     AND old_fragmentation.record_id = old_record.id
     AND old_record.title = new_record.title
       ;
""")
    

        
fragdb_partition_epilog = """

The `mmpdb partition` command takes one or more fragdb fragment files
created by `mmpdb fragment`, repartitions them by constant fragment,
and generates a new set of fragdb files that can be index in parallel.

The repartitioning places all fragmentations with the same constant
part into the same output fragdb file. It also includes a copy of all
of the structure records and error records in each output file. (The
structure records are needed for 1-cut hydrogen matching to work.)

The partitioning attempts to distribute the fragmentations evenly by
weight across the partition files. The weight is based on the number
of occurrences of each constant fragment. If there are `N` occurrences
then the weight is `N*(N-1)/2+1`.

There are two partitioning schemes. The `--num-files` / `-n` scheme
distributes the fragments using a first-fit bin packing algorithm
across at most `n` bins. The `--max-weight` scheme uses a first-fit
packing algorithm which distributes the constants across the available
bins such that no bin exceeds the given weight, with one exception; if
a constant's weight exceeds the maximum weight then the constant is
placed into a bin. If no bins are available then a new bin is added.

The constants file is a line-oriented format with one constant
fragment record per line. Each line contains white-space separated
columns. The first columns in the fragment SMILES. The second (used if
`--recount` is not specified) must be a positive integer for the
number of occurrences of the given fragment. If `--has-header` is
specified (the default) then the first line is treated as a header
lines and skipped. Use `--no-header` if there is no header.

Note: this is a superset of the output file generated by `mmpdb
fragdb_constants`.

If a constants file is specified then by default the constant SMILES
and and associated counts given in the file are used for the
partitioning. If `--recount` is specified then the counts are
determined from the specified fragdb files.

If a constants file is not specified then the constant SMILES and
associated counts are determined by analyzing the fragdb files.

The `--template` / `-t` string is used to generate the output
filenames, formatted using [Format String
Syntax](https://docs.python.org/3/library/string.html#formatstrings).
The available fields are:

\b
  {prefix} - the DATABASE path without the final extension
  {parent} - the DATABASE parent directory, or '.'
  {stem} - the DATABASE filename without the directory or final extension
  {sep} - the filesystem path seperator (eg, '/')
  {i} - an integer value 0 <= i < n

The value of `prefix`, `parent` and `stem` are determined from the
first fragdb DATABASE. For example, if the first DATABASE name is
'/abc/xyz.fragdb', on a macOS system, then the field values are:

\b
  {prefix} = '/abc/xyz'
  {parent} = '/abc'
  {stem} = 'xyz'
  {sep} = '/'

If only one DATABASE is specified then the default template is
`{prefix}-partition.{i:04}.fragdb` so the output filenames are based on
the input filename. For example:

\b
  % mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.fragdb -n 2 --dry-run
  i	#constants	weight	filename
  0	2	759121321	'ChEMBL_CYP3A4_hERG-partition.0000.fragdb'
  1	2	759116372	'ChEMBL_CYP3A4_hERG-partition.0001.fragdb'

If multiple DATABASES are specified, which may occur when the input
structures are fragmented in parallel into multiple fragdb files, then
the default template is `partition.{i:04}.fragdb` and the output
filenames are NOT based on the input filename. For example:

\b
  % mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb -n 2 --dry-run
  i	#constants	weight	filename
  0	2	759121321	'partition.0000.fragdb'
  1	2	759116372	'partition.0001.fragdb'

The `--dry-run` option has `fragdb_partition` print the partition
information to stdout but does not actually partition the file. This
lets you get an idea of how many partitions it will generate. (When
working with large data sets, it's best to generate and re-use a
constants file than have mmpdb determine the counts each time.)

Use `--tid` and `--num-tasks` to partition the files in parallel. The
idea is that `num_tasks` instances of `fragdb_partition` will be
submitted to a job scheduler and instance `tid` will partition only
the subsets `i` where `i % num_tasks == tid`.

"""


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
    help = "Maximum number of files to generate",
    )

@click.option(
    "--max-weight",
    type = positive_int(),
    default = None,
    callback = set_max_weight,
    help = (
        "Maximum weight per file (weight = N*(N-1)/2+1 where N "
        "is the number of occurrences of the constant"
        ),
    )


@click.option(
    "--constants",
    "-c",
    "constants_file",
    type = GzipFile("r"),
    help = "Only export fragmentations containing the constants specified in the named file",
    )

@click.option(
    "--recount",
    is_flag = True,
    help = "Ignore the counts in the constants file and instead compute them from fragment database(s)",
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
    type = template_type(),
    help = "Template for the output filenames",
    )
@click.option(
    "--tid",
    "--task-id",
    "task_id",
    type = nonnegative_int(),
    default = None,
    help = "Task id to use for parallelized partitioning (0 <= task id < num_tasks)",
    )
@click.option(
    "--num-tasks",
    type = positive_int(),
    default = None,
    help = "Number of tasks involved in parallelized partitioning",
    )
@click.option(
    "--dry-run",
    is_flag = True,
    default = False,
    help = "Describe the partitions which will be generated but don't generate them",
    )

@add_multiple_databases_parameters()
@click.pass_obj
def fragdb_partition(
        reporter,
        databases_options,
        num_files,
        max_weight,
        template,
        constants_file,
        recount,
        has_header,
        dry_run,
        task_id,
        num_tasks,
        ):
    """Partition fragments based on common constants

    DATABASE - one or more fragdb fragments database filenames
    """

    databases = databases_options.databases
    if not databases:
        raise click.BadArgumentUsage("must specify at least one fragment database")

    if task_id is None:
        if num_tasks is None:
            task_id = 0
            num_tasks = 1
        else:
            raise click.BadArgumentUsage("Must specify --task-id because --num-tasks was specified")
    else:
        if num_tasks is None:
            raise click.BadArgumentUsage("Must specify --num-tasks because --task-id was specified")
        if task_id >= num_tasks:
            raise click.BadArgumentUsage("--task-id must be less than --num-tasks")
    
    num_databases = len(databases)

    if template is None:
        if num_databases == 1:
            template = "{prefix}-partition.{i:04}.fragdb"
        else:
            template = "partition.{i:04}.fragdb"
    

    # Get values used for template generation
    reference = databases[0]
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
            constant_counts = get_constant_counts_from_file(
                constants_file, has_header, not recount)
        if not constant_counts:
            reporter.report("No constants found in {constants_file.name!r}. Exiting.")
            return            

        if recount:
            constant_counts = get_specified_constant_counts(list(constant_counts), databases, reporter)
        else:
            validate_databases(databases)

    # Sort from the most common to the least so we can use a first-fit algorithm.
    # (Decreasing count, increasing constant)
    constant_counts = list(constant_counts.items())
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
        

    if dry_run:
        click.echo(f"i\t#constants\tweight\tfilename")
        
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

        if not (subset_i % num_tasks == task_id):
            continue

        if dry_run:
            click.echo(f"{subset_i}\t{len(constant_subsets)}\t{total_weight}\t{output_filename!r}")
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
