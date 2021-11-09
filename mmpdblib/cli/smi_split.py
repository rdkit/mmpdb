import os
import click

from .click_utils import (
    command,
    template_type,
    positive_int,
    )

from .smi_utils import add_input_options

@command(name="smi_split")
@add_input_options

@click.option(
    "--num-files",
    "-n",
    default = None,
    type = positive_int(),
    help = "number of output SMILES files to generate",
    )

@click.option(
    "--num-records",
    default = None,
    type = positive_int(),
    help = "maximum number of SMILES to save to each file",
    )

@click.option(
    "--template",
    "-t",
    default = "{prefix}.{i:04}.smi",
    type = template_type(),
    show_default = True,
    help = "template for the output filenames",
    )

@click.argument(
    "smiles_filename",
    required = False,
    default = None,
    metavar = "FILE",
    )
@click.pass_obj
def smi_split(
        reporter,
        input_options,
        num_files,
        num_records,
        template,
        smiles_filename,
        ):
    """split the SMILES file 'FILE' into smaller files"""
    import pathlib
    from .. import fileio

    n = (num_files is not None) + (num_records is not None)
    if n == 2:
        # mutual exclusion is click is so tedious
        click.fail("Cannot specify both --num-files and --num-records")
    elif n == 0:
        num_files = 10


    # Get values used for template generation
    smiles_path = pathlib.Path("stdin.smi" if smiles_filename is None else smiles_filename)
    smiles_filename = str(smiles_path)
    smiles_parent = smiles_path.parent
    smiles_stem = smiles_path.stem
    smiles_prefix = str(smiles_path.parent / smiles_path.stem)

    # Read the SMILES
    try:
        with input_options.read_smiles_file(smiles_filename) as reader:
            entries = list(reader)

    except fileio.FileFormatError as err:
        die(f"Cannot parse input file: {err}")

    if not entries:
        reporter.report("No SMILES records found. Exiting.")
        return

    # Figure out how to partition
    if num_files is None:
        # use the number of records to get the number of files
        num_files = (len(entries) + num_records-1) // num_records
    else:
        # Can't have more files than entries
        num_files = min(num_files, len(entries))

    num_records_per_file = len(entries) // num_files

    i = -1
    for i in range(num_files):
        start = i * num_records_per_file
        end = start + num_records_per_file
        
        subset = entries[start:end]


        output_filename = template.format(
            parent = smiles_parent,
            stem = smiles_stem,
            prefix = smiles_prefix,
            sep = os.sep,
            i = i,
            )

        try:
            writer = fileio.open_output(output_filename, format_hint=None)
        except IOError as err:
            die(f"Cannot open output SMILES file: {err}")

        with writer:
            writer.writelines(f"{smiles}\t{id}\n" for (smiles, id) in subset)
    
    i += 1

    reporter.report(f"Created {num_files} SMILES files containing {len(entries)} SMILES records.")
    
