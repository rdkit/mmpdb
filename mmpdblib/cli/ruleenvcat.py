import click

from .click_utils import (
    command,
    GzipFile,
    add_single_database_parameters,
    open_dataset_from_options_or_exit,
    )

@command()
@click.option(
    "--pairs / --no-pairs",
    default = False,
    help = "With --pairs, include pairs in the output",
    )
@click.option(
    "--output",
    "-o",
    "outfile",
    default = "-",
    type = GzipFile("w"),
    help = "Write the rules to the named file (default is stdout)",
    )
@add_single_database_parameters()
@click.pass_obj
def ruleenvcat(
        reporter,
        database_options,
        pairs,
        outfile,
        ):
    "Show the rules in an mmpdb file"

    dataset = open_dataset_from_options_or_exit(database_options, quiet=True)
    rule_env_c = dataset.get_cursor()
    pair_c = dataset.get_cursor()

    
    outfile.write(f"from_smiles\tto_smiles\tradius\tSMARTS\tpseudoSMILES\n")
    for rule_env in dataset.iter_rule_environments(rule_env_c):
        outfile.write(
            f"{rule_env.from_smiles}\t{rule_env.to_smiles}\t{rule_env.radius}\t"
            f"{rule_env.smarts}\t{rule_env.pseudosmiles}\n"
            )
        if pairs:
            for pair in rule_env.iter_pairs(pair_c):
                outfile.write(f"\t{pair.compound1_id}\t{pair.compound2_id}\t{pair.constant_id}\n")
            
            
            
