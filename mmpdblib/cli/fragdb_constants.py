import math
import click
import collections
import contextlib

from .click_utils import (
    command,
    die,
    GzipFile,
    add_multiple_databases_parameters,
    open_fragdb_from_options_or_exit,
    positive_int,
    nonnegative_int,
    frequency_type,
)

COUNT_FRAGMENTATIONS_SQL = "SELECT COUNT(*) FROM fragmentation"
def _get_num_fragmentations(c):
    (num_fragmentations,) = next(c.execute(COUNT_FRAGMENTATIONS_SQL))
    return num_fragmentations
    

class SingleDatabase:
    def __init__(self, database):
        self.database = database
        self.db = self.c = None

    def __enter__(self):
        self.db = open_fragdb_from_options_or_exit(self.database)
        self.c = self.db.cursor()
        return self

    def __exit__(self, *args):
        self.c.close()
        self.db.close()
        self.db = self.c = None
        
    def get_num_fragmentations(self, reporter):
        return _get_num_fragmentations(self.c)

    def iter_constants(self, min_constant_num_heavies, min_count, max_count, reporter):
        query = """
  SELECT constant_smiles, n
    FROM (
      SELECT constant_smiles, count(*) AS n
        FROM fragmentation
       WHERE constant_num_heavies >= ?
    GROUP BY constant_smiles
    )
   WHERE ? <= n AND n <= ?
ORDER BY n DESC, constant_smiles
"""
        args = (min_constant_num_heavies, min_count, max_count)

        self.c.execute(query, args)
        return self.c
    

class MultipleDatabases:
    def __init__(self, databases):
        self.databases = databases

    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        pass
    
    def get_num_fragmentations(self, reporter):
        n = 0
        num_databases = len(self.databases)
        for database_i, database in enumerate(self.databases, 1):
            reporter.update(f"Analyzing {database!r} (#{database_i}/{num_databases})")
            try:
                with contextlib.closing(open_fragdb_from_options_or_exit(database)) as db:
                    with contextlib.closing(db.cursor()) as c:
                        n += _get_num_fragmentations(c)
            finally:
                reporter.update("")
        return n

    def iter_constants(self, min_constant_num_heavies, min_count, max_count,
                           reporter):
        # We need to load all of the constant SMILES counts
        query = """
   SELECT constant_smiles, count(*)
     FROM fragmentation
    WHERE constant_num_heavies >= ?
 GROUP BY constant_smiles
"""
        constant_counts = collections.Counter()
        num_databases = len(self.databases)
        for database_i, database in enumerate(self.databases, 1):
            reporter.update(f"Selecting constants from {database!r} (#{database_i}/{num_databases})")
            try:
                with contextlib.closing(open_fragdb_from_options_or_exit(database)) as db:
                    with contextlib.closing(db.cursor()) as c:
                        c.execute(query, (min_constant_num_heavies,))
                        for constant_smiles, n in c:
                            constant_counts[constant_smiles] += n
            finally:
                reporter.update("")
                        
        for constant_smiles, n in constant_counts.most_common():
            if min_count <= n <= max_count:
                yield constant_smiles, n

def open_frag_dbs(databases_options):
    databases = databases_options.databases
    if not databases:
        raise click.BadArgumentUsage("must specify at least one fragment database")
    if len(databases) == 1:
        return SingleDatabase(databases[0])
    else:
        return MultipleDatabases(databases)

                
fragdb_constants_epilog = """

By default this lists the constants in one or more fragdb files,
ordered from most common to least. It is meant as a way to reduce the
number of constants used during indexing, and to partition the
data set for parallel indexing.

Use `--min-count` and `--max-count` to set the minimum or maximum
number of occurences. (If a constant appears twice in the same record
then its occurence count is 2.)

Use `--min-frequency` and `--max-frequency` to express the minimum and
maximum occurences as a fraction of the total number of
occurences. (NOTE: this does not seem useful and will likely be
removed unless people say it's important.)

Use `--min-heavies-total-const-frag` to set a lower bound on the number
of heavies in each constant.

Use `--min-heavies-per-const-frag` to set a lower bound on the number
of heavies in the smallest fragment in each constant.

Use `--limit` to limite the output to the first `K` constants.

By default the constants are written to stdout. Use `--output` to
write the constants to a named file. If the filename ends with `.gz`
then the output is gzip compressed.

The output is formatted in two tab-separated columns as in the
following example:

\b
```
  % mmpdb fragdb_constants example.fragdb --limit 3
  constant	N
  *C	1010
  *C.*C	849
   *C.*O	662
```

The first column contains the fragment SMILES and the second contains
the count. The first line is a header with column named "constant" and
"N". Use `--no-header` to omit the header in the output.

"""


@command(
    epilog = fragdb_constants_epilog,
    name = "fragdb_constants",
    )
@click.option(
    "--min-count",
    type=nonnegative_int(),
)
@click.option(
    "--max-count",
    type=nonnegative_int(),
)
@click.option(
    "--min-frequency",
    "--min-freq",
    type=frequency_type(),
)
@click.option(
    "--max-frequency",
    "--max-freq",
    type=frequency_type(),
)
@click.option(
    "--min-heavies-per-const-frag",
    type=nonnegative_int(),
    help="Lower bound on the number of heavies in the smallest fragment in the constant part",
)
@click.option(
    "--min-heavies-total-const-frag",
    type=nonnegative_int(),
    default=0,
    help="Lower bound on the number of heavies in the constant part",
)
@click.option(
    "--limit",
    metavar="K",
    type=positive_int(),
    help="Limit the output to the 'K' most common constants",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    default = "-",
    type = GzipFile("w"),
    help = "Write the result to the named file (default: stdout)",
    )
@click.option(
    "--header / --no-header",
    default = True,
    help = "The default, --header, includes the header in output",
    )
@add_multiple_databases_parameters()
@click.pass_obj
def fragdb_constants(
    reporter,
    databases_options,
    min_count,
    max_count,
    min_frequency,
    max_frequency,
    min_heavies_per_const_frag,
    min_heavies_total_const_frag,
    limit,
    output_file,
    header,
):
    """List constants fragdb DATABASEs and their frequencies"""
    from ..index_algorithm import get_num_heavies

    with open_frag_dbs(databases_options) as frag_dbs:
        num_fragmentations = frag_dbs.get_num_fragmentations(reporter=reporter)

        const_frag_filter = min_heavies_per_const_frag is not None and min_heavies_per_const_frag > 0

        if min_count is None:
            min_count = 0

        if max_count is None:
            max_count = num_fragmentations

        # minimum frequency
        if min_frequency is not None:
            min_freq_count = int(math.ceil(min_frequency * num_fragmentations))
            min_count = max(min_count, min_freq_count)

        # maximum frequency
        if max_frequency is not None:
            max_freq_count = int(math.floor(max_frequency * num_fragmentations))
            max_count = min(max_count, max_freq_count)

        assert isinstance(min_count, int)
        assert isinstance(max_count, int)

        c = frag_dbs.iter_constants(min_heavies_total_const_frag, min_count, max_count,
                                        reporter=reporter)
        
        if limit is None:
            # 4611686018427387904 constants ought to be good enough for anyone
            limit = 2**63

        i = 0
        for constant_smiles, n in c:
            # Don't write the header until we have the first output line.
            # This makes the status reports easier to read as they are
            # not placed between the header and the constant lines.
            if header:
                output_file.write(f"constant\tN\n")
                # Only write the header once.
                header = False
            
            # Can't put the --limit in the SQL because of
            # possible additional filtering by
            # --min-heavies-per-constant-frag 
            # after the SQL query
            if i >= limit:
                break

            if const_frag_filter:
                terms = constant_smiles.split(".")
                if any(get_num_heavies(term) < min_heavies_per_const_frag for term in terms):
                    continue

            i += 1
            output_file.write(f"{constant_smiles}\t{n}\n")
    
    output_file.close()
