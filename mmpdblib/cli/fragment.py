"Implement the 'fragment' command"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2019-2021, Andrew Dalke Scientific, AB
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import click


from .. import smarts_aliases
from .click_utils import (
    command,
    die,
    name_to_command_line,
    ordered_group,
    positive_int,
)
from . import fragment_click
from . import smi_utils


# Make it so that ^C works in the main thread
def init_worker():
    import signal

    signal.signal(signal.SIGINT, signal.SIG_IGN)


def create_pool(num_jobs):
    import multiprocessing
    from .. import fragment_records

    if num_jobs > 1:
        pool = multiprocessing.Pool(num_jobs, init_worker)
    else:
        pool = fragment_records.SingleProcessPool()
    return pool


############ The fragment cli


fragment_epilog = """
Fragment molecules in a SMILES file by breaking on 'cut bonds', as
matched by --cut-smarts or the R-group SMILES of --cut-rgroup or
--cut-rgroup-file. Cut up to --num-cuts bonds. Don't fragment
molecules with more than --max-rotatable-bonds bonds or --max-heavies
heavy atoms. Don't create multiple cuts if the fragments in the
constant part have less than --min-heavies-per-const-frag atoms.  See
'mmpdb rgroup2smarts' for details about cutting with R-group SMILES.

The input structures come from a SMILES file. By default the fields
are whitespace delimited, where the first column contains the SMILES
and the second contains the identifier. Use --delimiter to change
the delimiter type, and --has-header to skip the first line. See
"mmpdb help-smiles-format" for more details.

The input SMILES strings are pre-processed to remove salts before
fragmenting. The default uses the default RDKit SmilesRemover. Use
--salt-remover to specify alternative rules, or use the special
name '<none>' to not remove any salts.

By default the fragmentation method uses 4 threads, which gives a
nearly 4-fold speedup. Use --num-jobs change the number of threads.

It can take a while to generate fragments. Suppose you want to update
the compound set every week, where only a few records are added or
changed. Most of the fragments will be the same between those two data
sets. What you can do is specify the old fragdb file as a --cache so
the fragmentation method can re-use the old fragmentation, assuming
the structure hasn't changed for a given record.

If you do not specify the `--output` filename then the default is
based on the SMILES filename with the "smi" or "smi.gz" extension
replace with "fragdb". If there is no input SMILES filename because
the data is from stdin then the default filename is `input.mmpdb`.

Examples:

1) Fragment the SMILES file to produce the fragdb file
`CHEMBL_thrombin_Ki_IC50.fragdb` (the output name is based on the
input SMILES filename):

\b
  % mmpdb fragment CHEMBL_thrombin_Ki_IC50.smi

2) Do the same, but with an explicit output filename:

\b
  % mmpdb fragment CHEMBL_thrombin_Ki_IC50.smi \\
      -o CHEMBL_thrombin_Ki_IC50.fragdb

3) Read from a gzip-compressed tab-delimited SMILES file. Use 8
threads to fragment the structures. Save the results to
dataset.fragdb:

\b
  % mmpdb fragment --delimiter tab dataset.smi.gz --num-jobs 8 \\
      -o dataset.fragdb

4) Fragment the SMILES in 'dataset.smi.gz'. Reuse fragment information
from the cache file 'old_dataset.fragdb' if possible, instead of
computing the fragments from scratch each time. Save the results to
'new_dataset.fragdb'.

\b
  % mmpdb fragment --cache old_dataset.fragdb dataset.smi.gz \\
      -o new_dataset.fragdb

\b
""" + smarts_aliases.get_epilog(
    "--cut-smarts", smarts_aliases.cut_smarts_aliases
)


def cannot_combine_with_fragment_options(ctx, cache):
    if cache is None:
        return
    used_names = ctx.meta[fragment_click.FRAGMENTATION_OPTION_NAMES]
    if not used_names:
        return
    names = sorted(name_to_command_line(name) for name in used_names)
    if len(names) == 1:
        raise click.UsageError(f"Cannot combine {names[0]} with --cache")
    else:
        *first, last = names
        first_str = ", ".join(first)
        raise click.UsageError(f"Cannot combine {first_str} or {last} with --cache")


@command(epilog=fragment_epilog)
@fragment_click.add_fragment_options
@click.option(
    "--cache",
    metavar="FRAGDB",
    help="Get fragment parameters and previous fragment information the FRAGDB file",
)
@click.option(
    "--num-jobs",
    "-j",
    metavar="N",
    type=positive_int(),
    default=4,
    help="Number of jobs to process in parallel (default: 4)",
)
@smi_utils.add_input_options
@click.option(
    "--output",
    "-o",
    metavar="FILENAME",
    help="Save the fragment data to FILENAME (default: based on the structure filename)",
)
@click.argument(
    "structure_filename",
    default=None,
    required=False,
    metavar="FILENAME",
    # help = "SMILES filename (default: read from stdin)",
)
@click.pass_context
def fragment(
    ctx,
    # from add_fragment_options
    fragment_options,
    # 'fragment'-specific arguments
    cache,
    num_jobs,
    # SMILES input options
    input_options,
    # output options
    output,
    # input
    structure_filename,
):
    """Fragment SMILES file structures on rotatable bonds

    FILENAME: SMILES file (default: read from stdin)

    The output is a 'fragdb' file containing fragmentations which can
    be used by `mmpdb index` or as cache for another `mmpdb fragment`.

    """
    from .. import (
        fragment_db,
        fragment_types,
        fragment_records,
        fileio,
    )

    config = ctx.obj
    cannot_combine_with_fragment_options(ctx, cache)

    output_filename = output
    if output_filename is None:
        if structure_filename is None:
            output_filename = "input.fragdb"
        else:
            output_filename = fileio.remove_suffixes(structure_filename) + ".fragdb"
        config.report(f"Using {output_filename!r} as the default --output file.")

    # Use a cache?
    cache_db = None
    if cache is not None:
        try:
            cache_db = fragment_db.open_fragdb(cache)
        except IOError as err:
            die(f"Cannot open cache: {err}")
        except ValueError as err:
            die(f"Problem loading cache: {err}")

        try:
            fragment_filter = cache_db.options.get_fragment_filter()
        except fragment_types.FragmentValueError as err:
            die(f"Error in cache option {err.name!r} ({err.value}!r) from {cache!r}: {err.reason}")
    else:
        try:
            fragment_filter = fragment_options.get_fragment_filter()
        except fragment_types.FragmentValueError as err:
            die("Error in command-line option %r (%r): %s" % (name_to_command_line(err.name), err.value, err.reason))

    pool = create_pool(num_jobs)

    try:
        try:
            with input_options.read_smiles_file(structure_filename) as reader:

                with fragment_db.open_fragment_writer(
                    output_filename,
                    options=fragment_filter.options,
                ) as writer:

                    records = fragment_records.make_fragment_records(
                        reader,
                        fragment_filter,
                        cache_db,
                        pool=pool,
                        reporter=config,
                    )
                    writer.write_records(records)

        except fileio.FileFormatError as err:
            die(f"Cannot parse input file: {err}")

        except UnicodeDecodeError as err:
            die(f"Error processing input file: {err} at {reader.location.where()}")

    except KeyboardInterrupt:
        config.update("Shutting down process pool")
        pool.terminate()
        pool.join()
        config.update("")
        raise SystemExit(-1)
    else:
        config.update("Closing process pool")
        pool.close()
        pool.join()
        config.update("")


####### "fragdb_utils" group


@ordered_group()
def fragdb_utils():
    pass


@fragdb_utils.command()
def fragdb_ls():
    pass


@fragdb_utils.command()
def fragdb_constant_stats():
    pass
