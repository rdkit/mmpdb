import math
import click

from .click_utils import (
    command,
    die,
    add_single_database_parameters,
    open_fragdb_from_options_or_exit,
    positive_int,
    nonnegative_int,
    frequency_type,
)

fragdb_constants_epilog = """
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
    help="lower bound on the number of heavies in the smallest fragment in the constant part",
)
@click.option(
    "--min-heavies-total-const-frag",
    type=nonnegative_int(),
    default=0,
    help="lower bound on the number of heavies in the constant part",
)
@click.option(
    "--limit",
    metavar="K",
    type=positive_int(),
    help="limit the output to the 'K' most common constants",
)
@add_single_database_parameters()
@click.pass_obj
def fragdb_constants(
    reporter,
    database_options,
    min_count,
    max_count,
    min_frequency,
    max_frequency,
    min_heavies_per_const_frag,
    min_heavies_total_const_frag,
    limit,
):
    """list constants in a fragdb DATABASE and their frequencies"""
    from ..index_algorithm import get_num_heavies

    fragdb = open_fragdb_from_options_or_exit(database_options)

    c = fragdb.cursor()

    (num_constants,) = next(c.execute("SELECT COUNT(*) FROM fragmentation"))

    # TODO: push this into the database?
    const_frag_filter = min_heavies_per_const_frag is not None and min_heavies_per_const_frag > 0

    if min_count is None:
        min_count = 0

    if max_count is None:
        max_count = num_constants

    # minimum frequency
    if min_frequency is not None:
        min_freq_count = int(math.ceil(min_frequency * num_constants))
        min_count = max(min_count, min_freq_count)

    # maximum frequency
    if max_frequency is not None:
        max_freq_count = int(math.floor(max_frequency * num_constants))
        max_count = min(max_count, max_freq_count)

    assert isinstance(min_count, int)
    assert isinstance(max_count, int)

    query = """
  SELECT constant_smiles, n
    FROM (
      SELECT constant_smiles, count(*) AS n
        FROM fragmentation
       WHERE constant_num_heavies >= ?
    GROUP BY constant_smiles
    )
   WHERE ? <= n AND n <= ?
ORDER BY n DESC
"""
    args = (min_heavies_total_const_frag, min_count, max_count)

    c.execute(query, args)

    if limit is None:
        # 4611686018427387904 constants ought to be good enough for anyone
        limit = 2**63

    click.echo(f"constant\tN")
    i = 0
    for constant_smiles, n in c:
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
        click.echo(f"{constant_smiles}\t{n}")
