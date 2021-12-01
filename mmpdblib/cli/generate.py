import click
from rdkit import Chem

from .. import environment
from ..index_algorithm import relabel
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
        
## class smiles_type(SMILES_MIXIN, click.ParamType):
##     name = "SMILES"

##     def convert(self, value, param, ctx):
##         # ignore None or molecule objects
##         if not isinstance(value, str):
##             return value
##         return self.convert_smiles(value, param, ctx)
    
class FragmentType(SMILES_MIXIN, click.ParamType):
    name = "SMILES"

    def __init__(self, add_atom_map):
        self.add_atom_map = add_atom_map
        
    def convert(self, value, param, ctx):
        # ignore None or molecule objects
        if not isinstance(value, str):
            return value

        if "*" not in value:
            self.fail("fragment SMILES must contain a '*'", param, ctx)
        if value.count("*") != 1:
            self.fail("fragent SMILES must contain only one '*'", param, ctx)
            
        mol = self.check_smiles(value, param, ctx)
        if mol.GetNumAtoms() < 2:
            self.fail("fragment SMILES must contain at least one heavy atom")

        if self.add_atom_map:
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atom.SetIntProp(ATOM_MAP_PROP, 1)
                    break

        return self.cansmi(mol)

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
    name = "hop",
    )
@click.option(
    "--constant",
    "constant_smiles",
    required = True,
    type = FragmentType(add_atom_map=True),
    help = "The constant SMILES",
    )

@click.option(
    "--query",
    "query_smiles",
    required = True,
    type = FragmentType(add_atom_map=False),
    help = "The query SMILES",
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
@add_single_database_parameters()
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="Explain each of the steps in the transformation process",
)

@click.pass_obj
def hop(
        reporter,
        query_smiles,
        constant_smiles,
        database_options,
        subqueries,
        radius,
        min_pairs,
        output_file,
        explain,
        ):
    reporter.set_explain(explain)
    
    dataset = open_dataset_from_options_or_exit(database_options, reporter.quiet)
    db = dataset.mmpa_db
    cursor = dataset.get_cursor()
    inner_cursor = dataset.get_cursor()

    # TODO: handle --subqueries
    query_smiles_list = [query_smiles]
    if subqueries:
        query_smiles_list.extend(get_subqueries(dataset, query_smiles, reporter))

    # From the constant part generate environment fingerprint
    reporter.explain(f"Using constant SMILES {constant_smiles} with radius {radius}.")
    unlabeled_constant_smiles = constant_smiles.replace("[*:1]", "*")
    centers = environment.find_centers(constant_smiles)
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

    # Find all rules with the given query environment and put them in a temporary table.
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
        # Get the rule_smiles.id for the query smiles
        labeled_query_smiles = relabel(query_smiles, [0])
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
            for cmpd1_id, cmdp1_smiles, cmpd2_id, cmpd2_smiles in inner_result:
                if swap:
                    cmpd1_id, cmdp1_smiles, cmpd2_id, cmpd2_smiles = cmpd2_id, cmpd2_smiles, cmpd1_id, cmdp1_smiles

                new_smiles, welded_mol = weld_fragments(unlabeled_constant_smiles, to_smiles)
                output_file.write(
                    f"{query_smiles}\t{new_smiles}\t{to_smiles}\t{num_pairs}\t"
                    f"{cmpd1_id}\t{cmdp1_smiles}\t{cmpd2_id}\t{cmpd2_smiles}\n"
                    )
                have_one = True
            assert have_one

        reporter.explain(f"Number of rules for {query_smiles}: {num_matching_rules}")
