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

from .. import fragment_types
from .. import fragment_records

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

FIELD_TITLES = {
    "start": (
        "start",
        "initial compound SMILES",
        ),
    "constant": (
        "constant",
        "constant fragment SMILES",
        ),
    "query": (
        "query",
        "query fragment SMILES",
        ),
    "radius": (
        "r",
        "environment radius",
        ),
    "smarts": (
        "smarts",
        "environment SMARTS",
        ),
    "pseudosmiles": (
        "pseudosmiles",
        "environment pseudo-SMILES",
        ),
    "generated": (
        "generated",
        "generated SMILES",
        ),
    "num_pairs": (
        "#pairs",
        "number of pairs from the corresponding rule",
        ),
    "lhs_id": (
        "lhs_id",
        "LHS id of a representative pair for the rule",
        ),
    "lhs_smiles": (
        "lhs_smiles",
        "LHS SMILES of a representative pair for the rule",
        ),
    "rhs_id": (
        "rhs_id",
        "RHS id of a representative pair for the rule",
        ),
    "rhs_smiles": (
        "rhs_smiles",
        "RHS SMILES of a representative pair for the rule",
        ),
    "rule_id": (
        "rule_id",
        "internal database rule id",
        ),
    "swapped": (
        "swapped",
        "direction of the transform in the rule",
        ),
    }

    
class ColumnFields(click.ParamType):
    name = "STR1,STR2,..."

    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value
        fields = value.split(",")
        for field in fields:
            if field not in FIELD_TITLES:
                known = ", ".join(repr(s) for s in FIELD_TITLES)
                self.fail(f"Unsupported column {field!r}. Available columns are {known}.")
        return fields
        

class ColumnTitles(click.ParamType):
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
    subqueries = []
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
            
        subqueries.append((-frag.constant_num_heavies, subquery))

    # Sort from largest to smallest
    subqueries.sort()
    subqueries = [subquery for (num_heavies, subquery) in subqueries]

    reporter.explain(f"Number of subqueries: {len(subqueries)}")
    if subqueries:
        reporter.explain(f"Subqueries are: {subqueries}")
    
    return subqueries


@command(
    name = "generate",
    )
@click.option(
    "--smiles",
    type = SmilesType(),
    help = "The full molecule",
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
    help = "If specified, also generate and include subqueries",
    )

@click.option(
    "--radius",
    default = 0, # XXX Is this right?
    type = radius_type(),
    help = "Fingerprint radius (default: 0)",
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
    default = "start,constant,query,radius,pseudosmiles,generated,num_pairs,lhs_id,lhs_smiles,rhs_id,rhs_id",
    show_default=True,
    help="a comma-separated list of output fields",
    )
@click.option(
    "--titles",
    type=ColumnTitles(),
    help="a comma-separated list of column titles (default based on the --fields)",
    )
    

@add_single_database_parameters()
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="Explain each of the steps in the transformation process",
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
        output_file,
        columns,
        titles,
        explain,
        ):
    """Generate 1-cut transforms"""
    reporter.set_explain(explain)

    num_options = (smiles is not None) + (query_smiles is not None) + (constant_smiles is not None)

    if (num_options == 0) or (num_options == 1 and smiles is None):
        raise click.UsageError("Must specify --smiles or two of --smiles, --constant, and --query")
    if num_options == 3:
        raise click.UsageError("Must specify --smiles or at most two of --smiles, --constant, and --query")

    if titles is None:
        titles = [FIELD_TITLES[column][0] for column in columns]
    else:
        if len(columns) != len(titles):
            raise click.UsageError(
                f"Mismatch between the number of --columns ({len(columns)}) and --titles ({len(titles)})"
                )
    title_line = "\t".join(titles) + "\n"
    output_format_str = "\t".join('{' + column + '}' for column in columns) + "\n"
    
    dataset = open_dataset_from_options_or_exit(database_options, reporter.quiet)
    cursor = dataset.get_cursor()
    inner_cursor = dataset.get_cursor()

    # Get fragmentation options from the database, limited to 1-cut
    options = dataset.get_fragment_options()
    options.num_cuts = 1
    fragment_filter = options.get_fragment_filter()

    if num_options == 1:
        assert smiles is not None
        queries = get_queries_from_smiles(smiles, fragment_filter, reporter)
        # Already enumerate all subqueries
        subqueries = False
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
        for constant_smiles, query_smiles_list in queries:
            assert len(query_smiles_list) == 1, (constant_smiles, query_smiles_list)
            query_smiles = query_smiles_list[0]
            query_smiles_list.extend(get_subqueries(dataset, query_smiles, reporter))

    output_file.write(title_line)
    for constant_smiles, query_smiles_list in queries:
        for result in generate_from_constant(
                        dataset = dataset,
                        cursor = cursor,
                        inner_cursor = inner_cursor,
                        constant_smiles = constant_smiles,
                        query_smiles_list = query_smiles_list,
                        radius = radius,
                        min_pairs = min_pairs,
                        reporter = reporter,
                        ):
            if 0:
                for k in result:
                    assert k in FIELD_TITLES, k
            output_line = output_format_str.format_map(result)
            output_file.write(output_line)
        

def generate_from_constant(
        dataset,
        cursor,
        inner_cursor,
        constant_smiles,
        query_smiles_list,
        radius,
        min_pairs,
        reporter,
        ):
    # I need to get the database.execute() so "?" is handled portably
    db = dataset.mmpa_db
    
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
        reporter.report("Constant SMILES environment fingerprint not present in the database.")
        return
    assert len(fpids) == 1, (env_fp.smarts, fpids)
    fpid = list(fpids)[0]

    # Find all rules with the given query environment
    result = db.execute("""
SELECT rule_id 
  FROM rule_environment
 WHERE environment_fingerprint_id = ?
""", (fpid,), cursor=cursor)

    rule_ids = [row[0] for row in result]
    if not rule_ids:
        reporter.warning(
            f"Found no rules using the SMILES environment fingerprint SMARTS {env_fp.smarts!r}.\n"
            "I don't think this should happen.")
        return

    reporter.explain(f"Number of matching environment rules: {len(rule_ids)}")


    for query_smiles in query_smiles_list:
        labeled_query_smiles = add_label_1(query_smiles)
        start_smiles, start_mol = weld_fragments(constant_smiles, labeled_query_smiles)
        
        # Get the rule_smiles.id for the query smiles
        labeled_query_smiles_id = dataset.get_rule_smiles_id(labeled_query_smiles, cursor=cursor)
        if labeled_query_smiles_id is None:
            reporter.explain(f"Query SMILES {labeled_query_smiles} is not a rule_smiles in the database.")
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
""", (labeled_query_smiles_id, labeled_query_smiles_id, fpid, min_pairs), cursor=cursor)

        num_matching_rules = 0
        for rule_id, from_smiles, to_smiles, rule_environment_id, num_pairs in result:
            num_matching_rules += 1
            #print(f"{query_smiles},{rule_id},{from_smiles},{to_smiles},{rule_environment_id},{num_pairs}")

            if from_smiles == labeled_query_smiles:
                swap = False
            elif to_smiles == labeled_query_smiles:
                from_smiles, to_smiles = to_smiles, from_smiles
                swap = True
            else:
                raise AssertionError(from_smiles, to_smiles, labeled_query_smiles)
            
            inner_result = db.execute("""
  SELECT cmpd1.public_id, cmpd1.clean_smiles, cmpd2.public_id, cmpd2.clean_smiles
    FROM pair,
         compound AS cmpd1,
         compound AS cmpd2
   WHERE pair.rule_environment_id = ?
     AND pair.compound1_id = cmpd1.id
     AND pair.compound2_id = cmpd2.id
ORDER BY cmpd1.clean_num_heavies * cmpd1.clean_num_heavies + cmpd2.clean_num_heavies * cmpd2.clean_num_heavies
 LIMIT 1
""", (rule_environment_id,), cursor = inner_cursor)
            have_one = False
            for lhs_id, lhs_smiles, rhs_id, rhs_smiles in inner_result:
                if swap:
                    lhs_id, lhs_smiles, rhs_id, rhs_smiles = rhs_id, rhs_smiles, lhs_id, lhs_smiles

                new_smiles, welded_mol = weld_fragments(constant_smiles, to_smiles)
                yield {
                    "start": start_smiles,
                    "constant": constant_smiles,
                    "query": query_smiles,
                    "radius": radius,
                    "smarts": env_fp.smarts,
                    "pseudosmiles": env_fp.pseudosmiles,
                    "generated": new_smiles,
                    "num_pairs": num_pairs,
                    "lhs_id": lhs_id,
                    "lhs_smiles": lhs_smiles,
                    "rhs_id": rhs_id,
                    "rhs_smiles": rhs_smiles,
                    "rule_id": rule_id,
                    "swapped": int(swap),
                    }
                have_one = True
            assert have_one

        reporter.explain(f"Number of rules for {query_smiles}: {num_matching_rules}")
