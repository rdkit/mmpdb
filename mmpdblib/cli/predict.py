"Implement the 'predict' command"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2019, Andrew Dalke Scientific, AB
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
    add_single_database_parameters,
    add_single_property,
    add_rule_selection_options,
    open_dataset_from_options_or_exit,
    )

# predict a property for a structure given the known property for
# another structure, using the MMPA database to identify the possible
# transformations and predicted property changes. The input to the
# tool will be the reference structure, the prediction structure, and
# the properties to predict. There will also be an option to report
# pairs with the same transformations and their properties.

predict_epilog = """

Predict the change in a specified '--property' going from the
'--reference' SMILES to the input '--smiles', using the matched
molecular pair 'DATABASE.' If several rules lead to the same input
'--smiles', use the one prioritized by '--rule-selection-cutoffs' and
'--score'.

By default the predicted delta value and standard deviation will be
written to stdout, in the format:
\b
  predicted delta: -4.91667 +/- 18.0477

If the optional reference property '--value' is specified then the
predicted output value will also be included, in the format: predicted
delta: -4.91667 predicted value: 7.08333 +/- 18.0477

The '--where' option takes a Python expression which describes which
transforms to consider. The '--rule-selection-cutoffs' and '--score'
options select which transform to use if multiple transforms produce
the same product. If the third-party APSW module is available then the
'--in-memory' option loads the SQLite database into memory before
doing the prediction. The '--explain' and '--times' options send
respectively debug and timing information to stderr. More details are
available from 'mmpdb transform --help'.

Use '--save-details' to save details about the possible transformations
to two tab-separated CSV files. The filenames are of the form
${PREFIX}_rules.txt and ${PREFIX}_pairs.txt, where ${PREFIX} is
specified by --prefix and by default is 'pred_details'.

The '${PREFIX}_rules.txt' file contains the following columns for each
transform, based on the rule, rule environment, and rule environment
statistics:
\b
  rule_id, rule_environment_id, radius, fingerprint, from_smiles,
  to_smiles, count, avg, std, kurtosis, skewness, min, q1, median, q3,
  max, paired_t, p_value

The '${PREFIX}_pairs.txt' file contains the following columns for
each of the pairs in each of the transforms:
\b
  from_smiles, to_smiles, radius, fingerprint, lhs_public_id,
  rhs_public_id, lhs_smiles, rhs_smiles, lhs_value, rhs_value, delta


Examples:

1) Predict the effect of substituting a sulfur in diphenyl ether if
the known melting point is 12C:

\b
  % mmpdb predict csd.mmpdb --smiles c1ccccc1Sc1ccccc1 \\
            --reference c1ccccc1Oc1ccccc1 --property MP --value 12.0
  predicted delta: +4.91667 predicted value: 16.9167 +/- 18.0477

2) Do the same calculation and save details about the transformations
to O_to_S_rules.txt and O_to_S_pairs.txt:

\b
  % mmpdb predict csd.mmpdb --smiles c1ccccc1Sc1ccccc1 \\
            --reference c1ccccc1Oc1ccccc1 --property MP --value 12.0 \\
            --save-details --prefix O_to_S

"""


# The following uses open_dataset_from_args_or_exit().
# In server use, if you have APSW, you may want to use:
#   dbutils.open_database(filename, copy_to_memory=True).get_dataset()
# This will copy the database into memory, instead of waiting for possible disk paging.
# (There may still be disk paging if you hit virtual memory limits.)


@command(epilog=predict_epilog)
@click.option(
    "--smiles",
    "-s",
    metavar="SMILES",
    required=True,
    help="The base structure to transform",
)
@click.option(
    "--reference",
    "reference_smiles",
    metavar="SMILES",
    required=True,
    help="The reference structure",
)
@add_single_database_parameters(add_in_memory=True)
@add_single_property
@add_rule_selection_options
@click.option(
    "--value",
    "-v",
    type=click.FLOAT,
    default=None,
    help="The property value for the reference",
)
@click.option(
    "--explain/--no-explain",
    help="Explain each of the steps in the transformation process",
)

## @click.option(
##     "--output",
##     "-o",
##     metavar = "FILENAME",
##     help = "save the output to FILENAME (default=stdout)",
##     )
@click.option(
    "--save-details",
    is_flag=True,
    default=False,
    help="Save information about the transformation pairs and statistics to two CSV files",
)
@click.option(
    "--prefix",
    metavar="STRING",
    default="pred_detail",
    help="Prefix to use for each CSV filename (default: 'pred_details')",
)
@click.option(
    "--times/--no-times",
    help="Report timing information for each step",
)
@click.pass_obj
def predict(
    reporter,
    database_options,
    rule_selection_options,
    property_name,
    explain,
    smiles,
    reference_smiles,
    value,
    save_details,
    prefix,
    times,
):
    """Predict the effect of a structural transformation

    DATABASE: an mmpdb database
    """
    import time
    from .. import (
        dbutils,
        analysis_algorithms,
    )

    reporter.set_explain(explain)
    # --smiles 'c1ccccc1C(=O)N(C)C' --reference 'c1ccccc1C(=O)NC' --property MP e.mmpdb --save-details
    start_time = time.time()

    dataset = open_dataset_from_options_or_exit(database_options, reporter.quiet)
    open_time = time.time()

    rule_selection_function = rule_selection_options.get_rule_selection_function()
    predict_tool = analysis_algorithms.get_predict_tool(dataset, rule_selection_function)

    if not predict_tool.is_available_property_name(property_name):
        property_names = sorted(predict_tool.get_property_names())
        if not property_names:
            die(f"--property {property_name!r} not available and no properties available")
        if len(property_names) == 1:
            die(f"--property {property_name!r} not available. Only {property_names[0]!r} is available.")

        available = ", ".join(map(repr, property_names))
        die(f"--property {property_name!r} not available. Available properties: {available}")

    reference_record = predict_tool.fragment_reference_smiles(reference_smiles)
    if reference_record.errmsg:
        die(f"Unable to fragment --reference {reference_smiles!r}: {reference_record.errmsg}")

    smiles_record = predict_tool.fragment_predict_smiles(smiles)
    if smiles_record.errmsg:
        die(f"Unable to fragment --smiles {smiles!r}: {reference_record.errmsg}")

    query_prep_time = time.time()

    explain = reporter.explain
    try:
        predict_result = predict_tool.predict(
            reference_record.fragmentations,
            smiles_record.fragmentations,
            property_name,
            explain=explain,
        )
    except analysis_algorithms.EvalError as err:
        sys.stderr.write("ERROR: %s\nExiting.\n" % (err,))
        raise SystemExit(1)

    predict_time = time.time()

    # Show the prediction
    property_rule = predict_result.property_rule
    if property_rule is None:
        print("No rules found.")
    else:
        delta = property_rule.avg
        std = property_rule.std

        # always include the "+" or "-" sign for a delta
        print_args = ["predicted delta: %+g" % (delta,)]
        if value is not None:
            print_args.append("predicted value: %g" % (value + delta,))
        if std is not None:
            print_args.append("+/- %g" % (std,))

        # Check if the rule exists in both directions
        if property_rule.is_bidirectional:
            print_args.append("(bidirectional)")

        print(*print_args)

    if save_details:
        try:
            property_rules_file = open(prefix + "_rules.txt", "w")
        except Exception as err:
            die("Unable to open output rules file: %s" % (err,))

        with property_rules_file:
            predict_result.write_property_rules(
                property_rules_file,
                field_names=(
                    "rule_environment_statistics_id",
                    "rule_id",
                    "rule_environment_id",
                    "radius",
                    "smarts",
                    "pseudosmiles",
                    "from_smiles",
                    "to_smiles",
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
                # column_aliases = {"rule_environment_id": "env_id"},
            )

        try:
            property_rule_pairs_file = open(prefix + "_pairs.txt", "w")
        except Exception as err:
            die("Unable to open output rules file: %s" % (err,))

        with property_rule_pairs_file:
            predict_result.write_property_rule_pairs(
                property_rule_pairs_file,
                field_names=(
                    "rule_environment_id",
                    "from_smiles",
                    "to_smiles",
                    "radius",
                    "smarts",
                    "pseudosmiles",
                    "lhs_public_id",
                    "rhs_public_id",
                    "lhs_smiles",
                    "rhs_smiles",
                    "lhs_value",
                    "rhs_value",
                    "delta",
                ),
                # column_aliases = {"from_smiles": "FROM_SMILES"},
            )
    output_time = time.time()

    if times:
        from .transform import get_time_delta_formatter

        sys.stderr.write("Elapsed time (in seconds):\n")
        format_dt = get_time_delta_formatter(output_time - start_time)

        sys.stderr.write("  open database: %s\n" % format_dt(open_time - start_time))
        sys.stderr.write("  prepare query: %s\n" % format_dt(query_prep_time - open_time))
        sys.stderr.write("        predict: %s\n" % format_dt(predict_time - query_prep_time))
        sys.stderr.write("   write output: %s\n" % format_dt(output_time - predict_time))
        sys.stderr.write("         TOTAL = %s\n" % format_dt(output_time - start_time))
