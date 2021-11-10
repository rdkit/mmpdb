import click

from .click_utils import (
    command,
    add_single_database_parameters,
    die,
    open_fragdb_from_options_or_exit,
    open_dataset_from_options_or_exit,
    )

@command(name="smicat")

@click.option(
    "--input-smiles",
    is_flag = True,
    default = False,
    help = "Use the input SMILES instead of the cleaned-up SMILES",
    )

@click.option(
    "--output",
    "-o",
    "output_file",
    default = "-",
    type = click.File("w"),
    help = "output filename (default is stdout)",
    )

@add_single_database_parameters()

def smicat(
        input_smiles,
        output_file,
        database_options,
        ):

    if database_options.database.endswith(".fragdb"):
        db = open_fragdb_from_options_or_exit(database_options)
        c = db.cursor()
        if input_smiles:
            c.execute("SELECT id, input_smiles FROM record UNION SELECT id, input_smiles FROM error_record")
        else:
            c.execute("SELECT id, normalized_smiles FROM record")
        iter_id_and_smiles = c
    else:
        dataset = open_dataset_from_options_or_exit(database_options)
        db = dataset.mmpa_db
        it = dataset.iter_compounds()
        if input_smiles:
            iter_id_and_smiles = ((compound.public_id, compound.input_smiles) for compound in it)
        else:
            iter_id_and_smiles = ((compound.public_id, compound.clean_smiles) for compound in it)

    with db:
        for id, smiles in iter_id_and_smiles:
            output_file.write(f"{smiles}\t{id}\n")

