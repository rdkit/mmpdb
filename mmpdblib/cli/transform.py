"Implement the 'transform' command"

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
import click
from .click_utils import (
    command,
    die,
    nonnegative_int,
    positive_int,
    radius_type,
    add_single_database_parameters,
    add_multiple_properties,
    add_rule_selection_options,
    get_property_names_or_error,
    open_dataset_from_options_or_exit,
)


class parse_smarts(click.ParamType):
    name = "SMARTS"

    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value

        from rdkit import Chem

        substructure_pat = Chem.MolFromSmarts(value)
        if substructure_pat is None:
            raise self.fail(f"Unable to parse SMARTS: {value!r}", param, ctx)
        return substructure_pat


# Helper function to make a new function which format time deltas
# so all the "."s are lined up and they use the minimal amount of
# left-padding.
def get_time_delta_formatter(max_dt):
    s = "%.1f" % (max_dt,)
    num_digits = len(s)
    fmt = "%" + str(num_digits) + ".1f"

    def format_dt(dt):
        return fmt % (dt,)

    return format_dt


transform_epilog = """

Apply transforms from an mmpdb database to an input structure.
Include possible property change statistics for the resulting
products.

Specify the input structure using --smiles. This will be fragmented
using the fragmentation parameters appropriate for the database. By
default all fragmentations will be considered. Use --min-variable-size
and --min-constant-size to set minimum heavy counts for the constant
and variable parts of the fragment.

By default the matching algorithm evaluates all radii around the local
environment of the constant's connection points. The scoring function
(see below) decides which radius is best. Use the --min-radius option
to require the environment match up to at least N bonds away.

Use --min-pairs to require that a transformation have at least N pairs
in the database. The default is 0, which allows all transformations.

For more complex filters, use the --where option. It takes a Python
expression which is allowed to use any of the following variable
names:

\b
  rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles,
  to_num_heavies, smirks, rule_environment_id, radius, fingerprint_id,
  fingerprint, rule_environment_statistics_id, count, avg, std,
  kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value,
  is_bidirectional

as well as the Python variables None, True, and False. The values come
from the transformation rule, the rule environment, and the rule
environment statistics.

The --where expression is evaluated independently for each
property. There is currently no way to specify a different expression
for different properties, or for the expression to know which property
is being evaluated.

Sometimes multiple transformations of an input structure lead to the
same final product. When this happens, the possible transformation are
put into bins based on their number of pairs. By default the first bin
contains transforms with at least 10 pairs, the second contains
transforms with at least 5 pairs, and the third contains transforms
with any pairs. The bin thresholds can be changed with the
--rule-selection-cutoffs option, which contains a list of
comma-separated integers. The default is equivalent to
'--rule-selection-cutoffs 10,5,0'.

A scoring function is used to decide which transformation to use from
a bin. The transformation with the largest score is selected. Use the
--score option to define an alternate scoring function as a Python
expression. The default is equivalant to:

  ((ninf if std is None else -std), radius, from_num_heavies, from_smiles)

This expression may contain any of the variables in the --where
option, as well as 'inf' and 'ninf', which can be used as equivalent
for positive and negative infinity in numerical comparisons. The
default expression selects the smallest standard deviation, breaks
ties using the largest radius, breaks ties with the most atoms
changed, then breaks ties arbitrarily by the substructure SMILES of
the left side of the transformation.

Note: if the --where or --score expressions start with a '-' then the
command-line parser may confuse it with a command-line option. In that
case, use a space as the first character in the expression, or enclose
the expression with parentheses.

Specify a SMARTS pattern with --substructure to limit the output to
products containing the specified substructure.

By default the output will contain predicted property changes for all
of the properties in the database. Use --property to select specific
options. Use it once for each property you want to evaluate. For
example, to get the results for both MW and MP use:
  --property MW --property MP

If you do not want property information in the output, specify
--no-properties.

The transformation output by default is sent to stdout. Use '--output'
to specify an output filename.

The output is a tab-delimited CSV file where the first line is a
header. The first column contains a sequential identifier, and the
second column contains the SMILES string for the product. Next come
the property columns. Each property gets 17 columns of output. The
column headers are prefixed with the property name followed by a "_"
and then the name of the data stored in the column. (There is
currently no way to limit the output to specific columns, other than
to change the code in do_transform.py.)

The "*_from_smiles" and "*_to_smiles" columns describe the
transformation.  Sometimes multiple transforms lead to the same
product and cause different properties to have different
transformations. The "*_radius" and "*_fingerprint" columns contain
the environment fingerprint radius and fingerprint string. The
remaining columns contain the database row id for the rule environment
record, and the statistics information for the given property.

The '--explain' option writes debug information to stderr. The
'--times' options reports timing information about the major stages to
stderr.

The transform code runs in a single thread by default. It has been
parallelized, though it is not highly scalable. Use '--jobs' to
specify the number of threads (really, "processes") to use.

Examples:

1) Generate all of the products of diphenyl ether using the MMP
transforms in the 'csd.mmpdb' database where there are are least 30
pairs. Also include the predicted effects on the 'MP' property. (Note:
the output is reformatted and trimmed for use as help text.)

\b
  % mmpdb transform csd.mmpdb --smiles 'c1ccccc1Oc1ccccc1' --min-pairs 15 -p MP
  ID                 SMILES MP_from_smiles       MP_to_smiles  MP_radius  \\
   1  COc1ccc(Oc2ccccc2)cc1  [*:1]c1ccccc1  [*:1]c1ccc(OC)cc1          0
   2             COc1ccccc1  [*:1]c1ccccc1             [*:1]C          0
   3   Cc1ccc(Oc2ccccc2)cc1  [*:1]c1ccccc1   [*:1]c1ccc(C)cc1          0
   4  Clc1ccc(Oc2ccccc2)cc1  [*:1]c1ccccc1  [*:1]c1ccc(Cl)cc1          0
   5              Oc1ccccc1  [*:1]c1ccccc1           [*:1][H]          0

\b
           MP_smarts  MP_pseudosmiles	MP_rule_environment_id	\\
  [#0;X1;H0;+0;!R:1]        [*:1](~*)                     1146
  [#0;X1;H0;+0;!R:1]        [*:1](~*)                      106
  [#0;X1;H0;+0;!R:1]        [*:1](~*)                      860
  [#0;X1;H0;+0;!R:1]        [*:1](~*)                     1199
  [#0;X1;H0;+0;!R:1]        [*:1](~*)                  1206787

\b
  MP_count   MP_avg  MP_std  MP_kurtosis  MP_skewness  MP_min   MP_q1  \\
        21	 16.333  32.273   2.3626         -1.2771      -22   -8.25
        28	-17.464  30.068	 -0.57255        -0.60008     -65  -41
        29	 5.5172  30.952	 -0.52187        -0.20078      45  -16.5
        20	16.3     33.304	 -0.41847         0.21611     -58   -6.5
        28	-6.4643  38.197	  0.51164         0.69676     -67  -32


\b
  MP_median  MP_q3  MP_max  MP_paired_t  MP_p_value
       10    32.75     116      -2.3193   0.03108
      -25.5   0.5       47       3.0734   0.0047958
        4    30         78      -0.9599   0.34532
       19    41         76      -2.1888   0.041303
      -10    15         99      -0.89551  0.37843

2) Require a standard deviation of no larger than 4.5 and give
priority to transformation with at least 20 pairs before following the
normal cutoffs. The --score here matches the default scoring function.

\b
  % mmpdb transform csd.mmpdb --smiles 'c1ccccc1Oc1ccccc1' \\
      --where 'std is not None and std < 4.5' \\
      --rule-selection-cutoffs '20,10,5,0' \\
      --score '((ninf if std is None else -std), radius, from_num_heavies, from_smiles)' \\
      --property MP
"""


@command(epilog=transform_epilog)
@add_single_database_parameters()
@click.option(
    "--smiles",
    "-s",
    required=True,
    help="The base structure to transform",
)
@click.option(
    "--min-variable-size",
    type=nonnegative_int(),
    metavar="N",
    default=0,
    help="Require at least N atoms in the variable fragment (default: 0)",
)
@click.option(
    "--max-variable-size",
    type=nonnegative_int(),
    default=9999,
    help="Allow at most N atoms in the variable fragment (default: 9999)",
)
@click.option(
    "--min-constant-size",
    type=nonnegative_int(),
    default=0,
    help="Require at least N atoms in the constant fragment (default: 0)",
)
@click.option(
    "--min-radius",
    "-r",
    type=radius_type(),
    default=0,
    help="Fingerprint radius (default: 0)",
)
@click.option(
    "--min-pairs",
    type=nonnegative_int(),
    default=0,
    help="Require at least N pairs in the transformation to report a product (default: 0)",
)
@click.option(
    "--substructure",
    "-S",
    "substructure_pat",
    type=parse_smarts(),
    help="Require the substructure pattern in the product",
)
@add_multiple_properties
@add_rule_selection_options
@click.option(
    "--jobs",
    "-j",
    "num_jobs",
    type=positive_int(),
    default=1,
    help="Number of jobs to run in parallel (default: 1)",
)
@click.option(
    "--explain",
    is_flag=True,
    default=False,
    help="Explain each of the steps in the transformation process",
)
@click.option(
    "--output",
    "-o",
    "output_filename",
    metavar="FILENAME",
    help="Save the output to FILENAME (default=stdout)",
)
@click.option(
    "--times",
    is_flag=True,
    default=False,
    help="Report timing information for each step",
)
@click.pass_obj
def transform(
    reporter,
    explain,
    min_radius,
    min_pairs,
    min_variable_size,
    max_variable_size,
    min_constant_size,
    database_options,
    property_names,
    no_properties,
    substructure_pat,
    rule_selection_options,
    smiles,
    num_jobs,
    output_filename,
    times,
    # **kwargs
):
    """Transform a structure

    DATABASE: a mmpdb database
    """
    import time
    from .. import (
        dbutils,
        analysis_algorithms,
        fileio,
    )

    assert max_variable_size > min_variable_size, "max-variable-size must be greater than min-variable-size"

    reporter.set_explain(explain)

    start_time = time.time()
    dataset = open_dataset_from_options_or_exit(database_options)
    open_time = time.time()

    property_names = get_property_names_or_error(dataset, property_names=property_names, no_properties=no_properties)

    if not property_names:
        include_empty = True
    else:
        include_empty = False  # should there be a --show-all option to enable this?

    # evaluate --where, --score, and --rule-selection-cutoffs.
    rule_selection_function = rule_selection_options.get_rule_selection_function()

    transform_tool = analysis_algorithms.get_transform_tool(dataset, rule_selection_function)
    transform_record = transform_tool.fragment_transform_smiles(smiles)

    if transform_record.errmsg:
        die(f"Unable to fragment --smiles {smiles!r}: {transform_record.errmsg}")

    transform_record = transform_tool.expand_variable_symmetry(transform_record)

    # Make sure I can open the output file before I start doing heavy work.
    try:
        outfile = fileio.open_output(output_filename, output_filename)
    except IOError as err:
        die(f"Cannot open --output file: {err}")

    query_prep_time = time.time()
    if num_jobs > 1:
        import multiprocessing

        pool = multiprocessing.Pool(processes=num_jobs)
    else:
        pool = None

    try:
        result = transform_tool.transform(
            transform_record.fragmentations,
            property_names,
            min_radius=min_radius,
            min_pairs=min_pairs,
            min_variable_size=min_variable_size,
            max_variable_size=max_variable_size,
            min_constant_size=min_constant_size,
            substructure_pat=substructure_pat,
            pool=pool,
            explain=reporter.explain,
        )
    except analysis_algorithms.EvalError as err:
        sys.stderr.write("ERROR: %s\nExiting.\n" % (err,))
        raise SystemExit(1)

    transform_time = time.time()

    with outfile:
        result.write_products(
            outfile,
            field_names=(
                #                "rule_environment_statistics_id",),
                "from_smiles",
                "to_smiles",
                "radius",
                "smarts",
                "pseudosmiles",
                "rule_environment_id",
                "count",
                "avg",
                "std",
                "kurtosis",
                "skewness",
                "min",
                "q1",
                "median",
                "q3",
                "max",
                "paired_t",
                "p_value",
            ),
            # column_aliases = {"from_smiles": "FROM"}, # use this to change the column name for a field
            include_empty=include_empty,
        )

    output_time = time.time()

    if times:
        sys.stderr.write("Elapsed time (in seconds):\n")
        format_dt = get_time_delta_formatter(output_time - start_time)
        sys.stderr.write("  open database: %s\n" % format_dt(open_time - start_time))
        sys.stderr.write("  prepare query: %s\n" % format_dt(query_prep_time - open_time))
        sys.stderr.write("      transform: %s\n" % format_dt(transform_time - query_prep_time))
        sys.stderr.write("   write output: %s\n" % format_dt(output_time - transform_time))
        sys.stderr.write("         TOTAL = %s\n" % format_dt(output_time - start_time))
