import click

from .click_utils import (
    command,
    add_single_database_parameters,
    die,
    GzipFile,
    open_fragdb_from_options_or_exit,
    open_dataset_from_options_or_exit,
    )

smicat_epilog = """

Each compound has two associated SMILES, the input SMILES used as
input to fragmentation, and the canonical SMILES string from RDKit
after input processing (typically desalting and structure
normalization). By default the output uses the processed SMILES. Use
`--input-smiles` to use the input SMILES string.

By default the output SMILES file is written to stdout. Use `--output`
to save the output to the named file.

Examples:

1) Write the cleaned-up structures as a SMILES file to stdout:

\b
  % mmpdb smicat csd.mmpdb

2) Save the structures to the file "original.smi", and use the input
SMILES instead of the de-salted SMILES:

\b
  % mmpdb smicat csd.mmpdb -o original.smi --input-smiles

"""

@command(
    name="smicat",
    epilog = smicat_epilog,
    )

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
    type = GzipFile("w"),
    help = "Output filename (default is stdout)",
    )

@add_single_database_parameters()

def smicat(
        input_smiles,
        output_file,
        database_options,
        ):
    """Write the mmpdb SMILES as a SMILES file"""

    if database_options.database.endswith(".fragdb"):
        db = open_fragdb_from_options_or_exit(database_options)
        c = db.cursor()
        if input_smiles:
            c.execute("SELECT title, input_smiles FROM record UNION SELECT title, input_smiles FROM error_record")
        else:
            c.execute("SELECT title, normalized_smiles FROM record")
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

