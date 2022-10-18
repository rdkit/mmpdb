import click

from .click_utils import (
    command,
    GzipFile,
    add_single_database_parameters,
    open_dataset_from_options_or_exit,
    )

@command()
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
def rulecat(
        reporter,
        database_options,
        outfile,
        ):
    "Show the rules in an mmpdb file"

    dataset = open_dataset_from_options_or_exit(database_options, quiet=True)
    rule_c = dataset.get_cursor()
    rule_env_c = dataset.get_cursor()

    outfile.write(f"id\tfrom_smiles\tto_smiles\n")
    for rule in dataset.iter_rules(rule_c):
        outfile.write(f"{rule.id}\t{rule.from_smiles}\t{rule.to_smiles}\n")
                    
    outfile.close()
