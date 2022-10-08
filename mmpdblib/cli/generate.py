import sys

import click
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
    return smiles.replace("*", "[*:1]")


class SMILES_MIXIN:
    def check_smiles(self, value, param, ctx):
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            self.fail("Cannot parse the SMILES", param, ctx)
        if get_num_frags(mol) != 1:
            self.fail("SMILES must not contain multiple fragments", param, ctx)

        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
            if atom.HasProp(ATOM_MAP_PROP):
                atom.ClearProp(ATOM_MAP_PROP)
            if atom.GetAtomicNum() == 0:
                if atom.GetFormalCharge() != 0:
                    self.fail("SMILES must not contain a '*' term with a charge", param, ctx)
                if atom.GetNumExplicitHs() != 0:
                    self.fail("SMILES must not contain a '*' term specified hydrogens", param, ctx)
                if len(list(atom.GetBonds())) != 1:
                    self.fail("SMILES '*' term must be a terminal atom", param, ctx)

        return mol

    def convert_smiles(self, value, param, ctx):
        mol = self.check_smiles(value, param, ctx)
        return self.cansmi(mol)

    def cansmi(self, mol):
        return Chem.MolToSmiles(mol)
        
class SmilesType(SMILES_MIXIN, click.ParamType):
    name = "SMILES"

    def convert(self, value, param, ctx):
        # ignore None or molecule objects
        if not isinstance(value, str):
            return value
        return self.convert_smiles(value, param, ctx)
    
class FragmentType(SMILES_MIXIN, click.ParamType):
    name = "SMILES"

    def convert(self, value, param, ctx):
        # ignore None or molecule objects
        if not isinstance(value, str):
            return value

        if "*" not in value:
            self.fail("fragment SMILES must contain a '*'", param, ctx)
        if value.count("*") != 1:
            self.fail("fragment SMILES must contain only one '*'", param, ctx)
            
        mol = self.check_smiles(value, param, ctx)
        if mol.GetNumAtoms() < 2:
            self.fail("fragment SMILES must contain at least one heavy atom")

        return self.cansmi(mol)

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
    return [(constant_smiles, [query_smiles])]

def get_queries_from_smiles_and_query(smiles, query_smiles, fragment_filter, reporter):
    for frag_term in get_queries_from_smiles(smiles, fragment_filter, reporter):
        frag_constant, frag_variables = frag_term
        if frag_variables[0] == query_smiles:
            return [frag_term]
    raise click.UsageError("--query SMILES not found in --smiles")

def get_queries_from_smiles_and_constant(smiles, constant_smiles, fragment_filter, reporter):
    for frag_term in get_queries_from_smiles(smiles, fragment_filter, reporter):
        frag_constant, frag_variables = frag_term
        if frag_constant == constant_smiles:
            return [frag_term]
    raise click.UsageError("--constant SMILES not found in --smiles")

###

def get_subqueries(dataset, query_smiles, reporter):
    # Get the subqueries SMILES by fragmenting the query fragment SMILES
    
    # Need to turn the "*" into something that can be processed.
    # Use argon as the placeholder.
    argon_smiles = query_smiles.replace("*", "[Ar]")

    # Need to use the fragmentation options from the database, limited to 1-cut
    options = dataset.get_fragment_options()
    options.num_cuts = 1
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

    # The fragmentations are symmetrical. I can limit my search to the constant_smiles.
    # There may be duplicate subqueries, so use a set to filter them out.
    subqueries = set()
    for frag in record.fragmentations:
        constant_smiles = frag.constant_smiles
        assert "[*]" not in constant_smiles
        if "[Ar]" in constant_smiles:
            if constant_smiles in ("*[Ar]", "[Ar]*"):
                # Ignore. Doesn't make sense.
                continue
            # Convert the [Ar] to an [H] and recanonicalize.
            mol = Chem.MolFromSmiles(constant_smiles.replace("[Ar]", "[H]"))
            assert mol is not None, constant_smiles
            subquery = Chem.MolToSmiles(mol)
        else:
            if frag.variable_smiles in ("*[Ar]", "[Ar]*"):
                # Isn't meaningful.
                continue
            subquery = constant_smiles
            
        subqueries.add((-frag.constant_num_heavies, subquery))

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
    type = FragmentType(),
    help = "The constant fragment SMILES",
    )

@click.option(
    "--query",
    "query_smiles",
    type = FragmentType(),
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
def generate(*args, **kwargs):
    prof_generate(*args, **kwargs)
    
def prof_generate(
        reporter,
        smiles,
        query_smiles,
        constant_smiles,
        database_options,
        subqueries,
        radius,
        min_pairs,
        output_file,
        columns,
        headers,
        include_header,
        num_jobs,
        chunksize,
        explain,
        ):
    """Apply 1-cut database transforms to a molecule"""
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

    # Figure out if I need to get a representative pair
    need_pair_colums = ("pair_from_id", "pair_from_smiles", "pair_to_id", "pair_to_smiles")
    if any((column in columns) for column in need_pair_colums):
        pair_cursor = dataset.get_cursor()
    else:
        pair_cursor = None

    # Get fragmentation options from the database, limited to 1-cut
    options = dataset.get_fragment_options()
    options.num_cuts = 1
    fragment_filter = options.get_fragment_filter()

    if num_options == 1:
        assert smiles is not None
        queries = get_queries_from_smiles(smiles, fragment_filter, reporter)
        # XXX Thinking about this
        ## subqueries = False
    else:
        if smiles is None:
            queries = get_queries_from_constant_and_query(
                constant_smiles, query_smiles)
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
        reporter.report("Single-threaded.")
    else:
        import multiprocessing
        pool = multiprocessing.Pool(processes=num_jobs)
        to_process = []
        reporter.report(f"Multiprocessing with {num_jobs} processors, chunksize={chunksize}.")

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
                
                    
            

def generate_unwelded_from_constant(
        dataset,
        cursor,
        pair_cursor,
        constant_smiles,
        from_smiles_list,
        radius,
        min_pairs,
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

    # Find all rules with the given query environment
    result = db.execute("""
SELECT COUNT(*)
  FROM rule_environment
 WHERE environment_fingerprint_id = ?
""", (fpid,), cursor=cursor)

    num_rule_ids = list(result)[0][0]
    if not num_rule_ids:
        reporter.warning(
            f"Found no rules using the SMILES environment fingerprint SMARTS {env_fp.smarts!r}.\n"
            "I don't think this should happen.")
        return

    reporter.explain(f"Number of matching environment rules: {num_rule_ids}")

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
SELECT t1.rule_id, from_smiles, to_smiles, rule_environment_id, n FROM
( SELECT rule.id AS rule_id,
         rule.from_smiles_id AS from_smiles_id,
         rule.to_smiles_id AS to_smiles_id,
         rule_environment.id AS rule_environment_id,
         count(*) as n
    FROM rule, rule_environment, pair
   WHERE (rule.from_smiles_id = ?
          OR rule.to_smiles_id = ?)
     AND environment_fingerprint_id = ?
     AND rule_environment.rule_id = rule.id
     AND pair.rule_environment_id = rule_environment.id
GROUP BY rule.id) t1
INNER JOIN
( SELECT rule.id AS rule_id,
         from_rule_smiles.smiles AS from_smiles,
         to_rule_smiles.smiles AS to_smiles
    FROM rule, rule_smiles AS from_rule_smiles, rule_smiles AS to_rule_smiles
   WHERE rule.from_smiles_id = from_rule_smiles.id 
     AND rule.to_smiles_id = to_rule_smiles.id) t2
ON t1.rule_id = t2.rule_id
   WHERE n >= ?
ORDER BY n DESC
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
                pair_result = db.execute("""
  SELECT cmpd1.public_id, cmpd1.clean_smiles, cmpd2.public_id, cmpd2.clean_smiles
    FROM pair,
         compound AS cmpd1,
         compound AS cmpd2
   WHERE pair.rule_environment_id = ?
     AND pair.compound1_id = cmpd1.id
     AND pair.compound2_id = cmpd2.id
ORDER BY cmpd1.clean_num_heavies * cmpd1.clean_num_heavies + cmpd2.clean_num_heavies * cmpd2.clean_num_heavies
 LIMIT 1
""", (rule_environment_id,), cursor = pair_cursor)
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
