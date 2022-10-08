"""Implement the 'index' command"""

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2021, Andrew Dalke Scientific AB
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

import sys
import os

import click

from .click_utils import (
    IntChoice,
    command,
    die,
    nonnegative_float,
    nonnegative_int,
    pop_known_args,
    positive_float,
    positive_int_or_none,
    set_click_attrs,
    open_fragdb_from_options_or_exit,
)

from .. import (
    config,
    index_types,
)


def get_psutil():
    try:
        import psutil
    except ImportError:
        psutil = None
    return psutil


#### Map click options to an index_types.IndexOptions


def add_index_options(command):
    DEFAULT_INDEX_OPTIONS = config.DEFAULT_INDEX_OPTIONS
    param_names = []

    def add_option(*args, **kwargs):
        # Keep track of the parameter names used
        param_names.append(args[0].lstrip("-").replace("-", "_"))

        click.option(*args, **kwargs)(command)

    add_option(
        "--min-variable-heavies",
        type=nonnegative_int(),
        default=DEFAULT_INDEX_OPTIONS.min_variable_heavies,
        help="Minimum number of non-hydrogen atoms in the variable fragment.",
    )

    add_option(
        "--max-variable-heavies",
        type=positive_int_or_none(),
        default=DEFAULT_INDEX_OPTIONS.max_variable_heavies,
        help=(
            "Maximum number of non-hydrogen atoms in the variable fragment "
            f"(default: {DEFAULT_INDEX_OPTIONS.max_variable_heavies}; for no maximum use 'none')"
        ),
    )

    add_option(
        "--min-variable-ratio",
        type=positive_float(),
        default=None,
        help="Minimum ratio of variable fragment heavies to heavies in the (cleaned) structure",
    )

    add_option(
        "--max-variable-ratio",
        type=nonnegative_float(),
        default=None,
        metavar="FLT",
        help="Maximum ratio of variable fragment heavies to heavies in the (cleaned) structure",
    )

    add_option(
        "--max-heavies-transf",
        type=nonnegative_int(),
        default=None,
        metavar="N",
        help="Maximum difference in the number of heavies transfered in a transformation",
    )

    add_option(
        "--max-frac-trans",
        type=nonnegative_float(),
        default=None,
        metavar="FLT",
        help="Maximum fraction of atoms taking part in a transformation",
    )

    add_option(
        "--min-radius",
        type=IntChoice(["0", "1", "2", "3", "4", "5"]),
        default=0,
        metavar="N",
        help="Minimum Environment Radius to be indexed in the MMPDB database",
    )

    add_option(
        "--max-radius",
        type=IntChoice(["0", "1", "2", "3", "4", "5"]),
        default=DEFAULT_INDEX_OPTIONS.max_radius,
        metavar="N",
        help="Maximum Environment Radius to be indexed in the MMPDB database",
    )

    assert DEFAULT_INDEX_OPTIONS.symmetric is False, "unsupported"
    add_option(
        "--symmetric",
        "-s",
        is_flag=True,
        default=DEFAULT_INDEX_OPTIONS.symmetric,
        help=(
            "Output symmetrically equivalent MMPs, i.e output both cmpd1,cmpd2, "
            "SMIRKS:A>>B and cmpd2,cmpd1, SMIRKS:B>>A"
        ),
    )

    assert DEFAULT_INDEX_OPTIONS.smallest_transformation_only is False, "unsupported"
    add_option(
        "--smallest-transformation-only",
        is_flag=True,
        default=DEFAULT_INDEX_OPTIONS.smallest_transformation_only,
        help="Ignore all transformations that can be reduced to smaller fragments",
    )

    def get_index_options_wrapper(**kwargs):
        # Fill in the defaults
        popped_kwargs = pop_known_args(param_names, kwargs, DEFAULT_INDEX_OPTIONS)

        # adjust for 'none'
        if popped_kwargs["max_variable_heavies"] == "none":
            popped_kwargs["max_variable_heavies"] = None

        check_validity(popped_kwargs)
        kwargs["index_options"] = index_types.IndexOptions(**popped_kwargs)

        return command(**kwargs)

    set_click_attrs(get_index_options_wrapper, command)

    return get_index_options_wrapper


def check_validity(kwargs):
    min_variable_heavies = kwargs["min_variable_heavies"]
    max_variable_heavies = kwargs["max_variable_heavies"]
    min_variable_ratio = kwargs["min_variable_ratio"]
    max_variable_ratio = kwargs["max_variable_ratio"]

    if min_variable_ratio is not None and max_variable_ratio is not None:
        if min_variable_ratio > max_variable_ratio:
            raise click.UsageError("--min-variable-ratio must not be larger than --max-variable-ratio")

    if min_variable_heavies is not None and max_variable_heavies is not None:
        if min_variable_heavies > max_variable_heavies:
            raise click.UsageError("--min-variable-heavies must not be larger than --max-variable-heavies")

    min_radius = kwargs["min_radius"]
    max_radius = kwargs["max_radius"]
    if min_radius > max_radius:
        raise click.UsageError("--min-radius must not be larger than --max-radius")

##### Get memory size

_process = None


def get_memory_use():
    global _process
    psutil = get_psutil()
    if psutil is None:
        return 0

    if _process is None:
        _process = psutil.Process(os.getpid())

    info = _process.memory_info()
    return info.rss  # or info.vms?


def human_memory(n):
    if n < 1024:
        return "%d B" % (n,)
    for unit, denom in (
        ("KB", 1024),
        ("MB", 1024 ** 2),
        ("GB", 1024 ** 3),
        ("PB", 1024 ** 4),
    ):
        f = n / denom
        if f < 10.0:
            return "%.2f %s" % (f, unit)
        if f < 100.0:
            return "%d %s" % (round(f, 1), unit)
        if f < 1000.0:
            return "%d %s" % (round(f, 0), unit)
    return ">1TB ?!?"


## for i in range(51):
##     print(2**i-1, human_memory(2**i-1))
##     print(2**i, human_memory(2**i))
##     print(2**i+1, human_memory(2**i+1))

################

index_epilog = """

Read a fragments file, match the molecular pairs, and save the results
to a mmpdb database file. Use "mmpdb help-analysis" for an explanation
about the terminology and process. A matched molecular pair connects
structure 1, made of variable part 1 (V1) and constant part C, with
structure 2, made of variable part 2 (V2) and constant part C, via the
transformation V1>>V2. The following uses |X| to mean the number of
heavy (non-isotopic hydrogen) atoms in X.

There are several ways to restrict which fragmentations or pairings
to allow. The --min-variable-heavies and --max-variable-heavies
options set limits to |V1| and |V2|. The --min-variable-ratio and
--max-variable-ratio set limits on |V1|/|V1+C| and |V2|/|V2+C|. The
--max-heavies-transf option sets a maximum bound on abs(|V1|-|V2|).

The --max-frac-trans places a limit on the number of atoms which take
part in the transformation, defined as |V+C(r)|/|V+C| where C(r) is
the number of atoms in the circular environment of the constant part
for a given radius r. C(0) is 0. The goal was to be able to exchange
transforms and environments but minimize the relative amount of
information revealed about the constant part.

The --max-radius option's default is set to 5 and can be changed to 
numbers less than 5, if you want to save DB space.

By default, mmpdb indexes all transformations per pair. To only index
one minimal transformation per pair, use
"--smallest-transformation-only". Note that this option is independent
of --max-radius, e.g. environment SMARTS/pseudoSMILES will still be
created for the smallest transformation.

The filter --max-variable-heavies is always used, with a default value
of 10. Specify "none" if you want no limit.

The transformation can be described as V1>>V2 or V2>>V1. By default
only one transformation is output; the one with the smallest value,
alphabetically. The --symmetric flag outputs both directions.

The optional --properties file contains physical property information
in a table format. Use "mmpdb help-property-format" for details. If a
property file is specified during the indexing stage then only those
compounds with at least one property are indexed. See "mmpdb
loadprops" for a way to specify properties after indexing.

By default the output will be saved in 'mmpdb' format using a filename
based on the input fragment filename and with the extension
'.mmpdb'. Use --output to specify an output filename. Use --out to
change the output format. By default the alternate formats will write
to stdout.

The 'mmpdb' format is based on a SQLite database. The 'mmpa' format
stores the same data in a text file with tab-separated fields. The
'csv' format is a tab-separated (not comma-separated!) table with the
columns:

\b
  SMILES1  SMILES2  id1  id2  V1>>V2  C

The 'csv' format does not include property information. The mmpa and csv
format also support gzip compression.

The --title specifies a string which goes into the database. The idea
is to store a label or short description to be displayed to the
user. If not given, it uses a short description based on the input
filename.

The --memory option writes a summary of memory use to stderr if the
'psutil' package is available. (If not, it prints a warning message
that it needs the package.) This was used during development to figure
out ways to reduce overall memory use.

The experimental --in-memory option loads the SQLite database into
memory before using it. This option requires the third-party APSW
SQLite adapater and will not work with Python's built-in SQLite
module. The transformation analysis does many scattered database
lookups. In some cases, like with a network file system, it can be
faster to load the entire database into memory where random access is
fast, rather than do a lot of disk seeks.


Examples:

1) Index 'csd.fragdb' and save the results to 'csd.mmpdb'. The
title will be "MMPs from 'csd.fragdb'".

\b
  % mmpdb index csd.fragdb

2) Index 'csd.fragdb', use properties from 'csd_MP.tsv' and limit the
pair matching to compounds listed in the TSV file, set the title to
'CSD MP', and save the results to 'csd_MP.mmpdb'.

\b
  % mmpdb index csd.fragdb --properties csd_MP.tsv --title "CSD MP" -o csd_MP.mmpdb

3) Limit the indexing to variable terms which have at least 12 heavy
atoms and where the size of the variable is no more than 40% of the
entire structure, and save the transformation in both A>>B and B>>A
(symmetric) forms:

\b
  % mmpdb index CHEMBL_thrombin_Ki_IC50.fragdb --symmetric \\
      --max-variable-ratio 0.4 --max-variable-heavies 12 \\
      --title "CHEMBL ratio 40%" --output CHEMBL_ratio_40.mmpdb
"""


@command(epilog=index_epilog)
@add_index_options
@click.option(
    "--properties",
    "properties_filename",
    metavar="FILENAME",
    help="File containing the identifiers to use and optional physical properties",
)
@click.option(
    "--output",
    "-o",
    "output_filename",
    metavar="FILENAME",
    help=(
        "Save the fragment data to FILENAME. "
        "Default for mmpdb is based on the fragment filename, "
        "otherwise stdout."
    ),
)
@click.option(
    "--out",
    "output_format",
    metavar="FORMAT",
    type=click.Choice([
        "csv", "csv.gz", "mmpa", "mmpa.gz", "mmpdb",
        "sql", "sql.gz", "sqlite", "sqlite.gz", "postgres", "postgres.gz",
        "csvd",
        ]),
    help=(
        "Output format. One of 'mmpdb' (default), 'csv', 'csv.gz', 'mmpa' or 'mmpa.gz'. "
        "If not present, guess from the filename, and default to 'mmpdb'"
    ),
)

@click.option(
    "--title",
    help="A short description of the dataset. If not given, base the title on the filename",
)
@click.option(
    "--memory",
    is_flag=True,
    default=False,
    help="Report a summary of the memory use",
)
@click.argument(
    "fragment_filename",
    metavar="FILENAME",
    type = click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
)
@click.pass_obj
def index(
    reporter,
    title,
    index_options,
    properties_filename,
    memory,
    output_format,
    output_filename,
    fragment_filename,
):
    """Index fragments and find matched molecular pairs

    FILENAME: the name of the fragdb file containing the fragments to index
    """
    from .. import (
        fragment_types,
        index_algorithm,
        properties_io,
    )

    if title is None:
        title = f"MMPs from {fragment_filename!r}"

    fragment_filter = index_options.get_fragment_filter()

    report_memory = False
    if memory:
        psutil = get_psutil()
        if psutil is None:
            sys.stderr.write("WARNING: Cannot report memory because the 'psutil' module is not available.\n")
        else:
            report_memory = True

    start_properties_memory = get_memory_use()
    if properties_filename is None:
        selected_ids = properties = None
    else:
        try:
            properties_file = open(properties_filename, "r", newline=None)
        except IOError as err:
            die(f"Cannot open --properties file: {err}")

        try:
            with properties_file:
                properties = properties_io.load_properties(properties_file, reporter)
        except ValueError as err:
            die(f"Problem reading --properties file {properties_file!r}: {err}")

        selected_ids = set(properties.get_ids())
    end_properties_memory = get_memory_use()

    if fragment_filename is None:
        # Not specified so use the default name.
        # XXX warning message?
        fragment_filename = "input.fragdb"

    if (output_format in (None, "mmpdb")) and output_filename is None:
        # Use the filename based on the fragments filename

        # replace the extension (if any) with ".mmpdb"
        output_filename = os.path.splitext(fragment_filename)[0] + ".mmpdb"
        reporter.warning(f"No --output filename specified. Saving to {output_filename!r}.")

    if (output_format == "csvd") and output_filename is None:
        # Use the directory name based on the fragments filename
        # replace the extension (if any) with ".csvd"
        output_filename = os.path.splitext(fragment_filename)[0] + ".csvd"
        reporter.warning(f"No --output directory specified. Saving to {output_filename!r}.")
        

    # reporter.report("Using fragment filters: %s" % (fragment_filter.get_args(),))

    start_fragment_index_memory = get_memory_use()
    fragment_reader = open_fragdb_from_options_or_exit(fragment_filename)

    with fragment_reader:
        with reporter.progress(fragment_reader, "Loaded fragment record") as report_fragment_reader:
            fragment_index = index_algorithm.load_fragment_index(report_fragment_reader, fragment_filter, selected_ids)

    start_mmp_memory = get_memory_use()
    environment_cache = index_algorithm.EnvironmentCache()

    pairs = index_algorithm.find_matched_molecular_pairs(
        fragment_index,
        fragment_reader,
        index_options,
        environment_cache,
        min_radius=index_options.min_radius,
        max_radius=index_options.max_radius,
        reporter=reporter,
    )

    pair_writer = index_algorithm.open_mmpa_writer(
        output_filename,
        format=output_format,
        title=title,
        fragment_options=fragment_reader.options,
        fragment_index=fragment_index,
        index_options=index_options,
        properties=properties,
        environment_cache=environment_cache,
        )
        
    with pair_writer:
        pair_writer.start()
        pair_writer.write_matched_molecule_pairs(pairs)
        end_mmp_memory = get_memory_use()
        pair_writer.end(reporter)
        end_memory = get_memory_use()

    if report_memory:
        # print(start_fragment_index_memory, start_mmp_memory, end_mmp_memory, end_memory)
        sys.stderr.write(
            "#pairs: %d Memory (RSS) total: %s properties: %s fragments: %s indexing: %s\n"
            % (
                pair_writer.num_pairs,
                human_memory(max(end_mmp_memory, end_memory)),
                human_memory(end_properties_memory - start_properties_memory),
                human_memory(start_mmp_memory - start_fragment_index_memory),
                human_memory(end_mmp_memory - start_mmp_memory),
            )
        )
