import click

from .click_utils import (
    command,
    add_multiple_databases_parameters,
    open_database_from_options_or_exit,
    )

@command()
@click.option(
    "--output",
    "-o",
    "outfile",
    default = "-",
    type = click.File("w"),
    )
@add_multiple_databases_parameters()
@click.pass_obj
def rulecat(
        reporter,
        databases_options,
        outfile,
        ):
    "Show the rules in an mmpdb file"
    from .. import dbutils

    for dbinfo, dataset in dbutils.iter_dbinfo_and_dataset(databases_options.databases, reporter):
        rule_c = dataset.get_cursor()
        rule_env_c = dataset.get_cursor()

        for rule in dataset.iter_rules(rule_c):
            outfile.write(f"{rule.from_smiles}>>{rule.to_smiles}\n")
                    
            
