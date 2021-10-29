import click

from .click_utils import command, IntChoice

# predict a property for a structure given the known property for
# another structure, using the MMPA database to identify the possible
# transformations and predicted property changes. The input to the
# tool will be the reference structure, the prediction structure, and
# the properties to predict. There will also be an option to report
# pairs with the same transformations and their properties.

predict_epilog="""

  --smiles SMILES, -s SMILES
                        the base structure to transform
  --reference SMILES    the reference structure
  --property NAME, -p NAME
                        property to use

  --value VALUE, -v VALUE
                        the property value for the reference

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
@add_single_dataset_arguments
@click.option(
    "--smiles",
    "-s",
    required = True,
    help = "the base structure to transform",
    )
@click.option(
    "--min-variable-size",
    metavar = "N",
    type = click.IntRange(0),
    default = 0,
    help = "require at least N atoms in the variable fragment (default: 0)",
    )

@click.option(
    "--max-variable-size",
    metavar = "N",
    type = click.IntRange(0),
    default = 9999,
    help = "allow at most N atoms in the variable fragment (default: 9999)",
    )

@click.option(
    "--min-constant-size",
    metavar = "N",
    type = nonnegative_int,
    default = 0,
    help = "require at least N atoms in the constant fragment (default: 0)",
    )

@click.option(
    "--min-radius",
    "-r",
    metavar = "R",
    type = IntChoice("0", "1", "2", "3", "4", "5"),
    default = 0,
    help = "fingerprint radius (default: 0)",
    )

@click.option(
    "--min-pairs",
    metavar = "N",
    type = click.IntRange(0),
    default = 0,
    help = "require at least N pairs in the transformation to report a product (default: 0)",
    )
    
p.add_argument(
    "--substructure",
    "-S",
    metavar = "SMARTS",
    help = "require the substructure pattern in the product",
    )
@add_multiple_properties_default_all
@add_rule_selection_arguments
@add_in_memory
@click.option(
    "--jobs",
    "-j",
    type = click.IntRange(1),
    default = 1,
    help = "number of jobs to run in parallel (default: 1)",
    )

@click.option(
    "--explain/--no-explain",
    help = "explain each of the steps in the transformation process",
    )

@click.option(
    "--output",
    "-o",
    metavar = "FILENAME",
    help = "save the output to FILENAME (default=stdout)",
    )
@click.option(
    "--times/--no-times",
    help = "report timing information for each step",
    )
def predict(**kwargs):
    "predict the effect of a structural transformation"
    # --smiles 'c1ccccc1C(=O)N(C)C' --reference 'c1ccccc1C(=O)NC' --property MP e.mmpdb --save-details
    start_time = time.time()
    
    dataset = dbutils.open_dataset_from_args_or_exit(args)
    open_time = time.time()
    
    rule_selection_function = analysis_algorithms.get_rule_selection_function_from_args(parser, args)
    predict_tool = analysis_algorithms.get_predict_tool(dataset, rule_selection_function)

    property_name = args.property
    if not predict_tool.is_available_property_name(property_name):
        property_names = sorted(predict_tool.get_property_names())
        if not property_names:
            parser.error("--property %r not available and no properties available"
                         % (property_name,))
        if len(property_names) == 1:
            parser.error("--property %r not available. Only %s is available."
                         % (property_name, property_names[0]))
        parser.error("--property %r not available. Available properties:"
                     % (property_name, ", ".join(property_names)))
            

    reference_record = predict_tool.fragment_reference_smiles(args.reference)
    if reference_record.errmsg:
        parser.error("Unable to fragment --reference %r: %s"
                     % (args.reference, reference_record.errmsg))

    smiles_record = predict_tool.fragment_predict_smiles(args.smiles)
    if smiles_record.errmsg:
        parser.error("Unable to fragment --smiles %r: %s"
                     % (args.smiles, reference_record.errmsg))
        
    query_prep_time = time.time()
    
    explain = command_support.get_explain(args.explain)
    try:
        predict_result = predict_tool.predict(
            reference_record.fragments, smiles_record.fragments, property_name,
            explain=explain)
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
        if args.value is not None:
            print_args.append("predicted value: %g" % (args.value + delta,))
        if std is not None:
            print_args.append("+/- %g" % (std,))

        # Check if the rule exists in both directions
        if property_rule.is_bidirectional:
            print_args.append("(bidirectional)")
    
        print(*print_args)

    if args.save_details:
        try:
            property_rules_file = open(args.prefix+"_rules.txt", "w")
        except Exception as err:
            parser.error("Unable to open output rules file: %s" % (err,))
            
        with property_rules_file:
            predict_result.write_property_rules(
                property_rules_file,
                field_names = (
                    "rule_environment_statistics_id",
                    "rule_id", "rule_environment_id", "radius", "fingerprint", "from_smiles", "to_smiles",
                    "count", "avg", "std", "kurtosis", "skewness",
                    "min", "q1", "median", "q3", "max", "paired_t", "p_value"),
                #column_aliases = {"rule_environment_id": "env_id"}, 
                    )
            
        try:
            property_rule_pairs_file = open(args.prefix+"_pairs.txt", "w")
        except Exception as err:
            parser.error("Unable to open output rules file: %s" % (err,))

        with property_rule_pairs_file:
            predict_result.write_property_rule_pairs(
                property_rule_pairs_file,
                field_names = (
                    "rule_environment_id",
                    "from_smiles", "to_smiles", "radius", "fingerprint",
                    "lhs_public_id", "rhs_public_id",
                    "lhs_smiles", "rhs_smiles",
                    "lhs_value", "rhs_value", "delta"),
                #column_aliases = {"from_smiles": "FROM_SMILES"},
                )
    output_time = time.time()

    if args.times:
        sys.stderr.write("Elapsed time (in seconds):\n")
        format_dt = get_time_delta_formatter(output_time - start_time)
        
        sys.stderr.write("  open database: %s\n" % format_dt(open_time - start_time))
        sys.stderr.write("  prepare query: %s\n" % format_dt(query_prep_time - open_time))
        sys.stderr.write("        predict: %s\n" % format_dt(predict_time - query_prep_time))
        sys.stderr.write("   write output: %s\n" % format_dt(output_time - predict_time))
        sys.stderr.write("         TOTAL = %s\n" % format_dt(output_time - start_time))
