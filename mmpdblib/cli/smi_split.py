import os
import click

from .click_utils import (
    command,
    template_type,
    positive_int,
    )

from .smi_utils import add_input_options

smi_split_epilog = """

Split a given SMILES file into several smaller SMILES files using one
of two available schemes. These new files can be used in a distributed
computing environment to fragment a dataset in parallel.

See `mmpdb help-smiles-format` for a description of how the
`--delimiter` and `--has-header` options affect parsing a SMILES
file. The output SMILES file has no header and uses tab-separated
columns.

The two split schemes read the input SMILES file and split the lines
(excluding any header) into the output files. The `--num-files` option
creates up to N output files in total, with roughly the same number of
records in each file. The `--num-records` option writes up to N
records to an output file, then switches to a new output file and
starts over.

If no scheme is specified then the default is the same as `--num-files
10`, which creates up to 10 output files.

The output filenames are based on a template pattern, which can be
changed with `--template`. The template may contain fields formatted
using Python's [Format String
Syntax](https://docs.python.org/3/library/string.html#formatstrings).

The available fields are:

\b
  {prefix} - the input SMILES filename without the final extension(s)
  {parent} - the parent directory of the input SMILES filename, or '.'
  {stem} - the SMILES filename without the directory or final extension
  {sep} - the filesystem path seperator (eg, '/')
  {i} - an integer value 0 <= i < n

The value of `prefix`, `parent` and `stem` are based on the input
SMILES filename. For example, if the filename is '/abc/xyz.smi', on
a macOS system, then the field values are:

\b
  {prefix} = '/abc/xyz'
  {parent} = '/abc'
  {stem} = 'xyz'
  {sep} = '/'

If the filename ends with '.gz' then that is also removed.

The default template is "{prefix}.{i:04}.smi". If the input filename
is `ChEMBL_CYP3A4_hERG.smi.gz` and there are 4 output files then the
output filenames are:

\b
```
ChEMBL_CYP3A4_hERG.0000.smi
ChEMBL_CYP3A4_hERG.0001.smi
ChEMBL_CYP3A4_hERG.0002.smi
ChEMBL_CYP3A4_hERG.0003.smi
```

If the output filename ends with ".gz", as with the template
"{prefix}.{i:04}.smi.gz", then the output files will be
gzip-compressed.

"""

@command(
    name="smi_split",
    epilog=smi_split_epilog)
@add_input_options

@click.option(
    "--num-files",
    "-n",
    default = None,
    type = positive_int(),
    help = "Number of output SMILES files to generate",
    )

@click.option(
    "--num-records",
    default = None,
    type = positive_int(),
    help = "Maximum number of SMILES to save to each file",
    )

@click.option(
    "--template",
    "-t",
    default = "{prefix}.{i:04}.smi",
    type = template_type(),
    show_default = True,
    help = "Template for the output filenames",
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
    """Split the SMILES file 'FILE' into smaller files"""
    import pathlib
    from .. import fileio

    n = (num_files is not None) + (num_records is not None)
    if n == 2:
        # mutual exclusion is click is so tedious
        click.fail("Cannot specify both --num-files and --num-records")
    elif n == 0:
        num_files = 10


    # Get values used for template generation
    s = "stdin.smi" if smiles_filename is None else smiles_filename
    # Trim '.gz' from the template
    if s.lower().endswith(".gz"):
        s = s[:-3]
    smi_path = pathlib.Path(s)
    smi_filename = str(smi_path)
    smi_parent = smi_path.parent
    smi_stem = smi_path.stem
    smi_prefix = str(smi_path.parent / smi_path.stem)

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
            parent = smi_parent,
            stem = smi_stem,
            prefix = smi_prefix,
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
    
