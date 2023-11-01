import sys

import click

class MolProcessingError(Exception):
    def __init__(self, error_message):
        super().__init__(error_message)
        self.error_message = error_message

from rdkit import Chem

from .. import environment
from ..analysis_algorithms import weld_fragments

from .click_utils import (
    command,
    die,
    radius_type,
    positive_int,
    GzipFile,
    add_single_database_parameters,
    open_dataset_from_options_or_exit,
    )

from .. import fragment_algorithm
from .. import fragment_records
from .. import fragment_types

ATOM_MAP_PROP = "molAtomMapNumber"

def get_num_frags(mol):
    return len(Chem.GetMolFrags(mol))

def add_label_1(smiles):
    return _add_label(smiles, 1)

def _add_label(smiles, i):
    left, mid, right = smiles.partition("*")
    if not mid:
        return smiles
    return f"{left}[*:{i}]{_add_label(right, i+1)}"

def parse_smiles_then_fragment(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise MolProcessingError("Cannot parse the SMILES")

    frags = Chem.GetMolFrags(mol)
    if not frags:
        raise MolProcessingError("No fragments found in the SMILES")
    
    return mol, frags

def check_no_wildcards(mol):
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        # We checked "*" earlier; this checks for "[#0]"
        raise MolProcessingError("SMILES must not contain a '*' term")


WHAT_OPTIONS = IS_MOLECULE, IS_CONSTANT, IS_VARIABLE = (101, 102, 103)

def sanitize_atom_properties(mol, frags, what):
    assert what in WHAT_OPTIONS
    
    num_frags = len(frags)
    if what == IS_MOLECULE:
        assert num_frags == 1, "wrong number of fragments"
        
    elif what == IS_CONSTANT:
        assert 1 <= num_frags <= 3, "wrong number of fragments"
        
    elif what == IS_VARIABLE:
        assert num_frags == 1, "wrong number of variable fragments"
        
    else:
        raise AssertionError("what?!")

    # Create a map from atom index to fragment index.
    # (Not needed when only one fragment, but simplifies the code.)

    atom_to_frag_idx = {}
    for frag_idx, atom_indices in enumerate(frags):
        for atom_idx in atom_indices:
            atom_to_frag_idx[atom_idx] = frag_idx
            
    fragments_with_wildcard = set()

    for atom in mol.GetAtoms():
        # Remove any isotope and atom map
        atom.SetIsotope(0)
        if atom.HasProp(ATOM_MAP_PROP):
            atom.ClearProp(ATOM_MAP_PROP)

        # Done processing non-wildcard atoms.
        if atom.GetAtomicNum() != 0:
            continue

        # Now we know we are processing a wildcard atom.

        if what == IS_MOLECULE:
            raise MolProcessingError("The SMILES may not include a '*' term")
        
        # If there is a "*", it must not have properties, and must be terminal
        
        if atom.GetFormalCharge() != 0:
            raise MolProcessingError("SMILES must not contain a '*' term with a charge")

        if atom.GetNumExplicitHs() != 0:
            raise MolProcessingError("SMILES must not contain a '*' term with a hydrogen count")

        if len(list(atom.GetBonds())) != 1:
            raise MolProcessingError("SMILES '*' term must be a terminal atom")
        
        # Ensure there is only one wildcard per fragment
        frag_idx = atom_to_frag_idx[atom.GetIdx()]
        if (what == IS_CONSTANT) and (frag_idx in fragments_with_wildcard):
            breakpoint()
            raise MolProcessingError("The constant SMILES may have only one '*' term on each fragment")

        fragments_with_wildcard.add(frag_idx)


    if what in (IS_CONSTANT, IS_VARIABLE) and len(fragments_with_wildcard) != num_frags:
        # Each fragment must have a wildcard
        if num_frags == 1:
            raise MolProcessingError("No '*' term found to indicate the attachment point")
        else:
            raise MolProcessingError("Each fragment must have a '*' term to indicate the attachment point")


# Used to prevent double processing of a click parameter.
class ProcessedSmiles:
    def __init__(self, canonical_smiles):
        self.canonical_smiles = canonical_smiles

class SmilesType(click.ParamType):
    name = "SMILES"

    def convert(self, value, param, ctx):
        # ignore None or ProcessedSmiles objects
        if not isinstance(value, str):
            return value

        # Do a quick check for "*"
        try:
            if "*" in value:
                raise MolProcessingError("SMILES must not contain a '*' term")

            mol, frags = parse_smiles_then_fragment(value)
            if len(frags) > 1:
                raise MolProcessingError("SMILES must contain only one fragment")

            sanitize_atom_properties(mol, frags, what = IS_MOLECULE)
            
        except MolProcessingError as err:
            self.fail(err.error_message, param, ctx)

        return ProcessedSmiles(Chem.MolToSmiles(mol))

    
class ConstantType(click.ParamType):
    name = "SMILES"

    def convert(self, value, param, ctx):
        # ignore None or ProcessedSmiles objects
        if not isinstance(value, str):
            return value

        try:
            mol, frags = parse_smiles_then_fragment(value)
            
            if not (1 <= len(frags) <= 3):
                raise MolProcessingError("Constant SMILES must contain only 1, 2, or 3 fragments")

            sanitize_atom_properties(mol, frags, what = IS_CONSTANT)

        except MolProcessingError as err:
            self.fail(err.error_message, param, ctx)
        
        return ProcessedSmiles(Chem.MolToSmiles(mol))
    
class QueryType(click.ParamType):
    name = "SMILES"

    def convert(self, value, param, ctx):
        # ignore None or ProcessedSmiles objects
        if not isinstance(value, str):
            return value

        try:
            mol, frags = parse_smiles_then_fragment(value)

            if len(frags) != 1:
                raise MolProcessingError("Query/variable SMILES must contain only one fragment")

            labeled_atoms = sanitize_atom_properties(mol, frags, what = IS_VARIABLE)
            
        except MolProcessingError as err:
            self.fail(err.error_message, param, ctx)
        
        # There may be a chirality loss when 1) there are two or three
        # wildcards, 2) the user indicates chirality, and 3) the
        # wildcards induce symmetry causing a center to disappear.

        # This is handled in mmpdb fragmentation by canonically
        # restoring chirality until the result matches the input.
        # That's not possible here since we only have the query.

        
        # I can think ways to restore it, eg, by replacing wildcards
        # with Sg, Lv and Db atoms, canonicalizing, restoring the
        # wildcards, recanonicalizing, and being careful about
        # mapping. I don't know if it will work though, and it sounds
        # like a lot of effort.

        # Up-enumeration should at least result in a close answer, so
        # the easy solution, and what I'll do, is to warn about a
        # difference in the number of chiral centers.

        def get_num_chiral_centers(s):
            num_double = s.count("@@")
            num_single = s.count("@") * 2*num_double
            return num_single + num_double

        new_smiles = Chem.MolToSmiles(mol)

        if new_smiles.count("*") > 1:
            old_count = get_num_chiral_centers(value)
            new_count = get_num_chiral_centers(new_smiles)
            if old_count != new_count:
                click.secho(
                    "Warning: Processing changed the number of query chiral centers from "
                    f"{old_count} to {new_count}.",
                    fg = "yellow",
                    err = True)

        return ProcessedSmiles(new_smiles)

COLUMN_DESCRIPTIONS = {
    "start": "initial compound SMILES (constant + from_smiles)",
    "constant": "constant fragment SMILES",
    "from_smiles": "SMILES of the fragment to modify (derived from the query)",
    "to_smiles": "SMILES of the modified fragment (from the database rule)",
    "r": "environment radius",
    "smarts": "environment SMARTS for the constant",
    "pseudosmiles": "environment pseudo-SMILES for the constant",
    "final": "generated SMILES (constant + to_smiles)",
    "heavies_diff": "change in the number of heavy atoms",
    "#pairs": "number of pairs from the environment rule",
    "pair_from_id": "id of the 'from' compound in a representative pair",
    "pair_from_smiles": "SMILES of the 'from' compound in a representative pair",
    "pair_to_id": "id of the 'to' compound in a representative pair",
    "pair_to_smiles": "SMILES of the 'to' compound in a representative pair",
    "rule_id": "internal database rule id",
    "env_rule_id": "internal database environment rule id",
    "swapped": "direction of the transform in the rule",
    }

COLUMN_EPILOG = "".join(f"\b\n* {column}: {description}" for (column, description) in COLUMN_DESCRIPTIONS.items())

DEFAULT_COLUMNS = (
        "start,constant,from_smiles,to_smiles,r,pseudosmiles,final,heavies_diff,"
        "#pairs,pair_from_id,pair_from_smiles,pair_to_id,pair_to_smiles"
        )
    
    
class ColumnFields(click.ParamType):
    name = "STR1,STR2,..."

    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value
        fields = value.split(",")
        for field in fields:
            if field not in COLUMN_DESCRIPTIONS:
                known = ", ".join(repr(s) for s in COLUMN_DESCRIPTIONS)
                self.fail(f"Unsupported column {field!r}. Available columns are {known}.")
        return fields
        

class ColumnHeaders(click.ParamType):
    name = "STR1,STR2,..."

    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value
        return value.split(",")

class SelectPairMethod(click.Choice):
    name = "method"
    def __init__(self):
        super().__init__(("first", "better", "quadratic", "min", "random"))

    def convert(self, value, param, ctx):
        value = super().convert(value, param, ctx)
        if value == "better":
            value = "quadratic"

        return value
    
def get_queries_from_smiles(smiles, fragment_filter, reporter):
    record = fragment_records.make_fragment_record_from_smiles(
        smiles,
        fragment_filter,
        reporter=reporter,
    )
    if record.errmsg:
        die(f"Cannot parse smiles: {record.errmsg}")

    queries = []
    for frag in record.fragmentations:
        queries.append((frag.constant_smiles, [frag.variable_smiles]))
    return queries

def get_queries_from_constant_and_query(constant_smiles, query_smiles):
    num_constant_attachments = constant_smiles.count("*")
    num_query_attachments = query_smiles.count("*")
    
    if num_constant_attachments != num_query_attachments:
        raise click.UsageError(
            "Mismatch between the number of attachment points in the --query "
            f"({num_query_attachments}) "
            "and the --constant "
            f"({num_constant_attachments})")

    return [(constant_smiles, [query_smiles])]

def get_queries_from_smiles_and_query(smiles, query_smiles, fragment_filter,
                                      reporter):
    reporter.explain(
        f"Fragmenting {smiles} to find the query {query_smiles}")
    
    for frag_term in get_queries_from_smiles(smiles, fragment_filter, reporter):
        frag_constant, frag_variables = frag_term
        if frag_variables[0] == query_smiles:
            reporter.explain(f"  => Found query. Constant is {frag_constant}")
            return [frag_term]
        reporter.explain(f"  Query does not match variable {frag_variables[0]}")
    raise click.UsageError("--query SMILES not found in --smiles")

def get_queries_from_smiles_and_constant(smiles, constant_smiles, fragment_filter,
                                        reporter):
    reporter.explain(
        f"Fragmenting {smiles} to find the constant {constant_smiles}")
    
    for frag_term in get_queries_from_smiles(smiles, fragment_filter, reporter):
        frag_constant, frag_variables = frag_term
        if frag_constant == constant_smiles:
            variable, = frag_variables
            reporter.explain(f"  => Found constant. Query is {variable}")
            return [frag_term]
        reporter.explain(f"  Constant does not match {frag_constant}")
    raise click.UsageError("--constant SMILES not found in --smiles")

###

def get_subqueries(dataset, query_smiles, reporter):
    # Get the subqueries SMILES by fragmenting the query fragment SMILES.
    # Replace the wildcards with argon atoms, fragment, and recover the subqueries.

    num_wildcards = query_smiles.count("*")
    assert num_wildcards > 0, query_smiles

    # Ignore cases where the constants are just the argons trimmed off.
    constants_to_ignore = {
        "*[Ar]", "*[Ar].*[Ar]", "*[Ar].*[Ar].*[Ar]",  # The canonical form RDKit generates.
        "[Ar]*", "[Ar].*[Ar]*", "[Ar]*.[Ar]*.[Ar]*",  # Makes me feel better to be complete.
        }
        
    # Turn the "*" into something that can be processed.
    # Use argon as the placeholder.
    argon_smiles = query_smiles.replace("*", "[Ar]")

    # Use the fragmentation options from the database.
    options = dataset.get_fragment_options()

    # Pointless to use more cuts than available wildcards.
    # (num_cuts sets the maximum number of cuts; cannot specify min_num_cuts.)
    options.num_cuts = min(options.num_cuts, num_wildcards)
    
    fragment_filter = options.get_fragment_filter()

    # Fragment
    record = fragment_records.make_fragment_record_from_smiles(
        argon_smiles,
        fragment_filter,
        reporter=reporter,
    )
    if not record.fragmentations:
        # No subqueries
        reporter.explain("The query SMILES cannot be fragmented using the database's fragmentation options.")
        return []

    # There may be duplicate subqueries, so use a set to filter them out.
    subqueries = set()
    
    for frag in record.fragmentations:
        if frag.num_cuts != num_wildcards:
            # Can't use as we need to match the number of input cuts.
            continue
                    
        variable_smiles = frag.variable_smiles
        assert "[*]" not in variable_smiles

        num_argons = variable_smiles.count("[Ar]")
        adjusted_variable_num_heavies = frag.variable_num_heavies - num_argons

        if num_argons > 0:
            # There is at least one argon in the variable part.
            
            # There must be at least one other, non-argon, heavy atom
            # otherwise we may end up with [H][H].
            if adjusted_variable_num_heavies == 0:
                continue

            # Convert the [Ar] atoms to [H] and recanonicalize.
            mol = Chem.MolFromSmiles(variable_smiles.replace("[Ar]", "[H]"))
            assert mol is not None, variable_smiles
            subquery = Chem.MolToSmiles(mol)

        else:
            # All of the argons have been cut off.

            # Since we are looking for *sub*queries, ignore the
            # trivial case where only the argons were cut.
            if frag.constant_smiles in constants_to_ignore:
                continue

            # Don't need to do more as the variable is already canonical.
            subquery = variable_smiles
            
        subqueries.add((-adjusted_variable_num_heavies, subquery))

    # Sort from largest to smallest and extract the SMILES
    subqueries = [subquery for (num_heavies, subquery) in sorted(subqueries)]

    reporter.explain(f"Number of subqueries: {len(subqueries)}")
    if subqueries:
        reporter.explain(f"Subqueries are: {subqueries}")
    
    return subqueries


def get_num_available_cpus():
    import os
    
    # The documentation for os.cpu_count() says:
    #   The number of usable CPUs can be obtained with
    #     ``len(os.sched_getaffinity(0))``
    # but that does not exist on my macOS/Python version

    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    
    # Fallback to the number of CPUs.
    #
    # This is what multiprocessing.Pool() does
    # when processes=None.
    
    n = os.cpu_count()

    if n is None:
        # number is indeterminable so be single threaded
        return 1
    
    return n

# The 'generate' command

generate_epilog = f"""

Use the DATABASE as a source of rules to generate new structures from
a given input structure. Only 1-cut rules are supported.

The core algorithm starts with a constant and a variable part (called
the "query"). The environment fingerprint for the constant, at radius
`-r` (default: 0) is found, and used to search the database for
matching rule environments where its corresponding rule has the query
SMILES on either side of the transform. If the environment has at
least `--min-pairs` pairs then use the other side of the tranform to
replace the query and generate a new structure.

The constant and variable parts can be specified directly, as
`--constant` and `--query`.

Alternatively, specify the `--smiles` and one of `--constant` or
`--query` and mmpdb will use the database's fragmentation options to
identify the other part.

Alternatively, if only `--smiles` is given then all of its constant
and variable parts will be used, again, fragmented according to the
database's fragmentation options.

The `--subqueries` option generates subfragments of the query
fragments and uses those as additional query fragments.

By default the output is written to stdout as tab-separated columns
with a header. Use `--no-header` to omit the header. Use `--output` to
write the output to a named file.

The following columns are available, though not all are included by
default:

{COLUMN_EPILOG}

Use `--columns`, containing a comma-separated list of column names, to
select a different set of columns or to change the output order. The
default is:

  {DEFAULT_COLUMNS}

By default the column header is the column name. Use `--headers` to
specify alternate names. This must be a comma-separated list of header
names, with the same length as `--columns`.

For example, the following generates columns with only the start
SMILES (with header "FROM_SMILES" and made from the subquery fragment
merged with the constant), and the final SMILES (with header
"TO_SMILES" and made from the SMILES on the other side of the matching
transform.

\b
```shell
  % mmpdb --quiet generate --smiles c1ccccc1CC1CC1 --radius 1 merged.mmpdb \\
      --subqueries --columns start,from_smiles,to_smiles,final \\
      --headers "START_SMILES,lhs_frag,rhs_frag,END_SMILES"
  START_SMILES	lhs_frag	rhs_frag	END_SMILES
  c1ccc(CC2CC2)cc1	[*:1]CC1CC1	[*:1][H]	c1ccccc1
  c1ccc(CC2CC2)cc1	[*:1]CC1CC1	[*:1]c1cc(C)nn1C	Cc1cc(-c2ccccc2)n(C)n1
  c1ccc(CC2CC2)cc1	[*:1]CC1CC1	[*:1]N1CCCCC1	c1ccc(N2CCCCC2)cc1
  c1ccc(C2CC2)cc1	[*:1]C1CC1	[*:1][H]	c1ccccc1
  c1ccc(C2CC2)cc1	[*:1]C1CC1	[*:1]OC	COc1ccccc1
     ... lines omitted ...
  CC1CC1	[*:1]C	[*:1]CCCC	CCCCC1CC1
  CC1CC1	[*:1]C	[*:1]COc1ccc(OC)cc1	COc1ccc(OCC2CC2)cc1
```

"""

@command(
    name = "generate",
    epilog = generate_epilog,
    )
@click.option(
    "--smiles",
    type = SmilesType(),
    help = "The full molecule to process",
    )

@click.option(
    "--constant",
    "constant_smiles",
    type = ConstantType(),
    help = "The constant fragment SMILES",
    )

@click.option(
    "--query",
    "--variable",
    "query_smiles",
    type = QueryType(),
    help = "The query/variable fragment SMILES",
    )

@click.option(
    "--subqueries / --no-subqueries",
    default = False,
    help = "If specified, also generate and include subfragments of the query fragment",
    )

@click.option(
    "--radius",
    default = 0, # XXX Is this right?
    type = radius_type(),
    help = "Fingerprint environment radius (default: 0)",
    )

@click.option(
    "--min-pairs",
    type = positive_int(),
    default = 1,
    help = "Only consider rules with at least N matched molecular pairs",
    )

@click.option(
    "--select-pair",
    "select_pair_method",
    type = SelectPairMethod(),
    default = "first",
    help = (
        "If 'first' (fastest), select a representative pair arbitrarily. "
        "If 'quadratic' or 'better', minimize sum of num_heavies**2. "
        "If 'min', use the minimum num_heavies for either side. "
        "If 'random', select one at random."
        ),
    )

@click.option(
    "--output",
    "-o",
    "output_file",
    default = "-",
    type = GzipFile("w"),
    )
@click.option(
    "--columns",
    type=ColumnFields(),
    default = DEFAULT_COLUMNS,
    help=f"A comma-separated list of output fields (see below for the default)",
    )
@click.option(
    "--headers",
    type=ColumnHeaders(),
    help="A comma-separated list of column headers (default uses --fields)",
    )
@click.option(
    # Not using --header since it's to easy to get mixed up with --headers
    "--no-header",
    "include_header",
    is_flag=True,
    default=True,
    help = "Use --no-header to exclude the column headers in the output",
    )
    
@click.option(
    "--num-jobs",
    "-j",
    default = 1,
    type = click.IntRange(0),
    help = "Number of processes to use when welding SMILES (0 means use all available CPUs)",
    )

@click.option(
    "--chunksize",
    default = 100,
    type = click.IntRange(1),
    help = "Number of SMILES to process in each multiprocessing work unit (default: 100)",
    )
    
@add_single_database_parameters(add_in_memory=True)
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="Explain the steps in the generation process",
)

@click.pass_obj
def generate(
        reporter,
        smiles,
        query_smiles,
        constant_smiles,
        database_options,
        subqueries,
        radius,
        min_pairs,
        select_pair_method,
        output_file,
        columns,
        headers,
        include_header,
        num_jobs,
        chunksize,
        explain,
        ):
    """Apply database transforms to a molecule"""
    reporter.set_explain(explain)

    num_options = (smiles is not None) + (query_smiles is not None) + (constant_smiles is not None)

    if (num_options == 0) or (num_options == 1 and smiles is None):
        raise click.UsageError("Must specify --smiles or two of --smiles, --constant, and --query")
    if num_options == 3:
        raise click.UsageError("Must specify --smiles or at most two of --smiles, --constant, and --query")

    if headers is None:
        headers = columns[:]
    else:
        if len(columns) != len(headers):
            raise click.UsageError(
                f"Mismatch between the number of --columns ({len(columns)}) and --headers ({len(headers)})"
                )
    header_line = "\t".join(headers) + "\n"
    output_format_str = "\t".join('{' + column + '}' for column in columns) + "\n"
    
    dataset = open_dataset_from_options_or_exit(database_options, reporter.quiet)
    cursor = dataset.get_cursor()
#    cursor.execute("PRAGMA mmap_size=2147418112")

    # Figure out if I need to get a representative pair
    need_pair_colums = ("pair_from_id", "pair_from_smiles", "pair_to_id", "pair_to_smiles")
    if any((column in columns) for column in need_pair_colums):
        pair_cursor = dataset.get_cursor()
    else:
        pair_cursor = None

    # Get fragmentation options from the database.
    options = dataset.get_fragment_options()
    fragment_filter = options.get_fragment_filter()

    # Unwrap the canonical SMILES
    if smiles is not None:
        smiles = smiles.canonical_smiles
    if query_smiles is not None:
        query_smiles = query_smiles.canonical_smiles
    if constant_smiles is not None:
        constant_smiles = constant_smiles.canonical_smiles
        
    if num_options == 1:
        assert smiles is not None
        queries = get_queries_from_smiles(smiles, fragment_filter, reporter)
    else:
        if smiles is None:
            queries = get_queries_from_constant_and_query(constant_smiles, query_smiles)
            
        elif constant_smiles is None:
            queries = get_queries_from_smiles_and_query(
                smiles, query_smiles, fragment_filter, reporter)
            
        elif query_smiles is None:
            queries = get_queries_from_smiles_and_constant(
                smiles, constant_smiles, fragment_filter, reporter)
            
        else:
            raise AssertionError

    if subqueries:
        for constant_smiles, from_smiles_list in queries:
            assert len(from_smiles_list) == 1, (constant_smiles, from_smiles_list)
            from_smiles = from_smiles_list[0]
            from_smiles_list.extend(get_subqueries(dataset, from_smiles, reporter))

    if include_header:
        output_file.write(header_line)

    if num_jobs == 0:
        num_jobs = get_num_available_cpus()
        
    if num_jobs == 1:
        import contextlib
        pool = contextlib.nullcontext()
        to_process = None
    else:
        import multiprocessing
        pool = multiprocessing.Pool(processes=num_jobs)
        to_process = []

    with pool:
        for constant_smiles, from_smiles_list in queries:
            for unwelded_result in generate_unwelded_from_constant(
                            dataset = dataset,
                            cursor = cursor,
                            pair_cursor = pair_cursor,
                            constant_smiles = constant_smiles,
                            from_smiles_list = from_smiles_list,
                            radius = radius,
                            min_pairs = min_pairs,
                            select_pair_method = select_pair_method,
                            reporter = reporter,
                            ):
                if to_process is None:
                    result = weld_unwelded_result(unwelded_result)
                    output_line = output_format_str.format_map(result)
                    output_file.write(output_line)
                else:
                    to_process.append(unwelded_result)

        if to_process is not None:
            reporter.report(f"Welding {len(to_process)} results.")
            for result in pool.imap_unordered(
                    weld_unwelded_result, to_process, chunksize=chunksize):
                output_line = output_format_str.format_map(result)
                output_file.write(output_line)
                

SELECT_PAIR_SQL = """
 SELECT cmpd1.public_id, cmpd1.clean_smiles, cmpd2.public_id, cmpd2.clean_smiles
   FROM pair,
        compound AS cmpd1,
        compound AS cmpd2
  WHERE pair.rule_environment_id = ?
    AND pair.compound1_id = cmpd1.id
    AND pair.compound2_id = cmpd2.id
<ORDER>
  LIMIT 1
"""
_select_pair_sql_table = {
    "first": SELECT_PAIR_SQL.replace(
        "<ORDER>\n",
        ""),
    "quadratic": SELECT_PAIR_SQL.replace(
        "<ORDER>",
        "ORDER BY cmpd1.clean_num_heavies * cmpd1.clean_num_heavies + cmpd2.clean_num_heavies * cmpd2.clean_num_heavies"),
    "min": SELECT_PAIR_SQL.replace(
        "<ORDER>",
        "ORDER BY MIN(cmpd1.clean_num_heavies, cmpd2.clean_num_heavies)"),
    "random": SELECT_PAIR_SQL.replace(
        "<ORDER>",
        "ORDER BY random()"),
    }
            

def generate_unwelded_from_constant(
        dataset,
        cursor,
        pair_cursor,
        constant_smiles,
        from_smiles_list,
        radius,
        min_pairs,
        select_pair_method,
        reporter,
        ):
    # I need to get the database.execute() so "?" is handled portably
    db = dataset.mmpa_db

    # Set some defaults
    pair_from_id = pair_from_smiles = pair_to_id = pair_to_smiles = None
    
    # From the constant part generate environment fingerprint
    reporter.explain(f"Using constant SMILES {constant_smiles} with radius {radius}.")
    labeled_constant_smiles = add_label_1(constant_smiles)
    centers = environment.find_centers(labeled_constant_smiles)
    env_fps = environment.compute_constant_environment_from_centers(centers, radius, radius)
    assert len(env_fps) == 1, (constant_smiles, centers, env_fps)
    env_fp = env_fps[0]
    reporter.explain(f"Environment SMARTS: {env_fp.smarts} pseudoSMILES: {env_fp.pseudosmiles}")
    
    # Get the environment fingerprint id
    fpids = dataset.get_smarts_ids([env_fp.smarts], cursor=cursor)
    if not fpids:
        reporter.explain("Constant SMILES environment fingerprint not present in the database.")
        return
    assert len(fpids) == 1, (env_fp.smarts, fpids)
    fpid = list(fpids)[0]

    ## Commented out since it takes a lot of time and doesn't seem useful.
    # (1.6 seconds on one benchmark, out of 4.9 seconds total.)
##     # Find all rules with the given query environment
##     result = db.execute("""
## SELECT COUNT(*)
##   FROM rule_environment
##  WHERE environment_fingerprint_id = ?
## """, (fpid,), cursor=cursor)

##     num_rule_ids = list(result)[0][0]
##     if not num_rule_ids:
##         reporter.warning(
##             f"Found no rules using the SMILES environment fingerprint SMARTS {env_fp.smarts!r}.\n"
##             "I don't think this should happen.")
##         return

##     reporter.explain(f"Number of matching environment rules: {num_rule_ids}")

    for from_smiles in from_smiles_list:
        labeled_from_smiles = add_label_1(from_smiles)
        start_smiles, start_mol = weld_fragments(constant_smiles, labeled_from_smiles)
        start_num_heavies = start_mol.GetNumHeavyAtoms()
        
        # Get the rule_smiles.id for the query smiles
        labeled_from_smiles_id = dataset.get_rule_smiles_id(labeled_from_smiles, cursor=cursor)
        if labeled_from_smiles_id is None:
            reporter.explain(f"Query SMILES {labeled_from_smiles} is not a rule_smiles in the database.")
            continue

        # Find rules with the given SMILES on either side and enough pairs
        result = db.execute("""
SELECT t1.rule_id,
       from_rule_smiles.smiles,
       to_rule_smiles.smiles,
       t1.rule_environment_id,
       t1.n
FROM
(
SELECT rule.id AS rule_id,
       rule.from_smiles_id as from_smiles_id,
       rule.to_smiles_id as to_smiles_id,
       rule_environment.id AS rule_environment_id,
       rule_environment.num_pairs AS n
    FROM
        rule,
        rule_environment

   WHERE (rule.from_smiles_id = ? OR
          rule.to_smiles_id = ?)
     AND environment_fingerprint_id = ?
     AND rule_environment.rule_id = rule.id
     AND rule_environment.num_pairs >= ?
) t1,
  -- Also include the from-/to- SMILES
  rule_smiles AS from_rule_smiles,
  rule_smiles AS to_rule_smiles
   WHERE t1.from_smiles_id = from_rule_smiles.id
     AND t1.to_smiles_id = to_rule_smiles.id

   ORDER BY n DESC, from_rule_smiles.smiles, to_rule_smiles.smiles
""", (labeled_from_smiles_id, labeled_from_smiles_id, fpid, min_pairs), cursor=cursor)

        num_matching_rules = 0
        for rule_id, from_smiles, to_smiles, rule_environment_id, num_pairs in result:
            num_matching_rules += 1
            #print(f"{from_smiles},{rule_id},{from_smiles},{to_smiles},{rule_environment_id},{num_pairs}")

            if from_smiles == labeled_from_smiles:
                swap = False
            elif to_smiles == labeled_from_smiles:
                from_smiles, to_smiles = to_smiles, from_smiles
                swap = True
            else:
                raise AssertionError(from_smiles, to_smiles, labeled_from_smiles)

            ## Welding takes most of the time so don't do the calculation here.
            # new_smiles, welded_mol = weld_fragments(constant_smiles, to_smiles)
            # final_num_heavies = fragment_algorithm.count_num_heavies(welded_mol)

            if pair_cursor is not None:
                # Pick a representative pair
                select_pair_sql = _select_pair_sql_table[select_pair_method]
                    
                pair_result = db.execute(
                    select_pair_sql,
                    (rule_environment_id,),
                    cursor = pair_cursor,
                    )
                
                have_one = False
                for pair_from_id, pair_from_smiles, pair_to_id, pair_to_smiles in pair_result:
                    have_one = True
                    if swap:
                        pair_from_id, pair_from_smiles, pair_to_id, pair_to_smiles = (
                            pair_to_id, pair_to_smiles, pair_from_id, pair_from_smiles
                            )
                assert have_one

            yield {
                "start": start_smiles,
                "constant": constant_smiles,
                "from_smiles": from_smiles,
                "to_smiles": to_smiles,
                "r": radius,
                "smarts": env_fp.smarts,
                "pseudosmiles": env_fp.pseudosmiles,
                #"final": new_smiles,
                #"heavies_diff": final_num_heavies - start_num_heavies,
                "start_num_heavies": start_num_heavies,
                "#pairs": num_pairs,
                "pair_from_id": pair_from_id,
                "pair_from_smiles": pair_from_smiles,
                "pair_to_id": pair_to_id,
                "pair_to_smiles": pair_to_smiles,
                "rule_id": rule_id,
                "swapped": int(swap),
                }

        reporter.explain(f"Number of rules for {from_smiles}: {num_matching_rules}")

def weld_unwelded_result(d):
    constant_smiles = d["constant"]
    to_smiles = d["to_smiles"]
    start_num_heavies = d.pop("start_num_heavies")
    
    new_smiles, welded_mol = weld_fragments(constant_smiles, to_smiles)
    
    final_num_heavies = welded_mol.GetNumHeavyAtoms()
    d["final"]= new_smiles
    d["heavies_diff"] = final_num_heavies - start_num_heavies
    return d
