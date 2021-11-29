import click
from rdkit import Chem

from .. import environment
from ..index_algorithm import relabel

"""
Input 
constant = [*]c1ccccc1
query = [*]C1CC1

2) Once we have a query environment fingerprint: search in the "environment fingerprint table" in the mmpdb database. 
If the query environment fingerprint is present, then get its index from the environment fingerprint table. If not present 
then search ends here. Let's assume that we found a query environment fingerprint in the database and its index is 100.

3) search [*]C1CC1 in the "rule_smiles table". If [*]C1CC1 is present in the rule_smiles table, get its index.
If it's not present then search ends here. Let's assume that we found [*]C1CC1 in rule_smiles table and its index is 5

4) Search the query environment using its index i.e. 100 in rule_envionment table and get all matching rules. Here we get rule ids
Lets assume we find four rules whose indexes are: 1, 2, 3, 4

5) For each of the four rule ids, search in the "rule table" and check if  [*]C1CC1 (use its index i.e 5) is present in either LHS or RHS side.
 Let's assume we found some rules where  [*]C1CC1 is present in either LHS or RHS. The output at this step might look like this:

envfp_id/ rule_id/LHS/RHS
100/1/5/10
100/3/5/12
100/4/20/5

5) Now we need to get smiles for new fragments from their indexes which are 10, 12 and 20. Look inside the rule_smiles table
and get smiles. Let's assume for 10 its [*]Cl, for 12 its [*]C1CCC1 and for 20 its [*]Br

6) Take a constant part from the query molecule i.e [*]c1ccccc1 and weld it with [*]Cl, [*]Br and [*]C1CCC1.
This will produce three compounds:
Clc1ccccc1
Brc1ccccc1
C1CC(C1)c2ccccc2

"""
from .click_utils import (
    command,
    die,
    radius_type,
    positive_int,
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


@command(
    name = "hop",
    )
@click.option(
    "--constant",
    "constant_smiles",
    required = True,
    type = FragmentType(add_atom_map=True),
    help = "the constant SMILES",
    )

@click.option(
    "--query",
    "query_smiles",
    required = True,
    type = FragmentType(add_atom_map=False),
    help = "the query SMILES",
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
    help = "fingerprint radius (default: 0)",
    )

@click.option(
    "--min-pairs",
    type = positive_int(),
    default = 1,
    help = "only consider rules with at least N matched molecular pairs",
    )

@click.option(
    "--output",
    "-o",
    "output_file",
    default = "-",
    type = click.File("w"),
    )
@add_single_database_parameters()
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="explain each of the steps in the transformation process",
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

    ## Get the variable SMILES by fragmenting the input SMILES and
    ## finding the specified constant.

    ## # Need to use the fragmentation options from the database, limited to 1-cut
    ## options = dataset.get_fragment_options()
    ## options.num_cuts = 1
    ## fragment_filter = options.get_fragment_filter()

    ## # Fragment
    ## record = fragment_records.make_fragment_record_from_smiles(
    ##     smiles,
    ##     fragment_filter,
    ##     reporter=reporter,
    ## )
    ## if not record.fragmentations:
    ##     die("The input SMILES cannot be fragmented using the database's fragmentation options.")

    ## # Find the matching fragmentation
    ## variable_smiles = None
    ## for frag in record.fragmentations:
    ##     if frag.constant_smiles == constant_smiles:
    ##         variable_smiles = frag.variable_smiles
    ##         break
    ## else:
    ##     lines = [
    ##         f"No fragmentation found using the constant SMILES {constant_smiles!r}.",
    ##         "The available constants are:",
    ##         ]
    ##     for frag in record.fragmentations:
    ##         lines.append(f"  {frag.constant_smiles}   with variable: {frag.variable_smiles}")
        
    ##     die(*lines)

    ## reporter.explain(f"Using variable {variable_smiles!r}.")

    # TODO: handle --subqueries
    query_smiles_list = [query_smiles]

    # From the constant part generate environment fingerprint
    reporter.explain(f"Using constant SMILES {constant_smiles} with radius {radius}.")
    centers = environment.find_centers(constant_smiles)
    env_fps = environment.compute_constant_environment_from_centers(centers, radius, radius)
    assert len(env_fps) == 1, (constant_smiles, centers, env_fps)
    env_fp = env_fps[0]
    reporter.explain(f"Environment SMARTS: {env_fp.smarts}")
    
    
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

    # Add the 

    for query_smiles in query_smiles_list:
        # Get the rule_smiles.id for the query smiles
        labeled_query_smiles = relabel(query_smiles, [0])
        labeled_query_smiles_id = dataset.get_rule_smiles_id(labeled_query_smiles, cursor=cursor)
        if labeled_query_smiles_id is None:
            reporter.explain(f"Query SMILES {labeled_query_smiles} is not a rule_smiles in the database.")
            continue

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
        for rule_id, from_smiles, to_smiles, rule_environment_id, num_pairs in result:
            #print(f"{query_smiles},{rule_id},{from_smiles},{to_smiles},{rule_environment_id},{num_pairs}")
            
            
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
                print(
                    f"{query_smiles}\t{from_smiles}\t{to_smiles}\t{num_pairs}\t"
                    f"{cmpd1_id}\t{cmdp1_smiles}\t{cmpd2_id}\t{cmpd2_smiles}"
                    )
                have_one = True
            assert have_one
