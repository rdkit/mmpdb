# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
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

from __future__ import print_function

import sys
import time
import multiprocessing

from rdkit import Chem

from . import command_support
from . import dbutils
from . import analysis_algorithms
from . import fileio

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

########################

def transform_command(parser, args):
    min_radius = args.min_radius
    assert min_radius in list("012345"), min_radius
    min_radius = int(min_radius)
    min_pairs = int(args.min_pairs)
    min_variable_size = args.min_variable_size
    max_variable_size = args.max_variable_size
    assert max_variable_size > min_variable_size, "max-variable-size must be greater than min-variable-size"
    min_constant_size = args.min_constant_size
    
    explain = command_support.get_explain(args.explain)

    start_time = time.time()
    dataset = dbutils.open_dataset_from_args_or_exit(args)
    open_time = time.time()
    
    property_names = command_support.get_property_names_or_error(parser, args, dataset)
    if not property_names:
        include_empty = True
    else:
        include_empty = False  # should there be a --show-all option to enable this?

    if args.substructure:
        substructure_pat = Chem.MolFromSmarts(args.substructure)
        if substructure_pat is None:
            parser.error("Cannot parse --substructure %r" % (args.substructure,))
    else:
        substructure_pat = None
        
    # evaluate --where, --score, and --rule-selection-cutoffs.
    rule_selection_function = analysis_algorithms.get_rule_selection_function_from_args(
        parser, args)
    
    transform_tool = analysis_algorithms.get_transform_tool(dataset, rule_selection_function)
    transform_record = transform_tool.fragment_transform_smiles(args.smiles)
    transform_record = transform_tool.expand_variable_symmetry(transform_record)

    if transform_record.errmsg:
        parser.error("Unable to fragment --smiles %r: %s"
                     % (args.smiles, transform_record.errmsg))

    # Make sure I can open the output file before I start doing heavy work.
    try:
        outfile = fileio.open_output(args.output, args.output)
    except IOError as err:
        parser.error("Cannot open --output file: %s" % (err,))

    query_prep_time = time.time()
    if args.jobs > 1:
        pool = multiprocessing.Pool(processes=args.jobs)
    else:
        pool = None
    try:
        result = transform_tool.transform(
            transform_record.fragments, property_names,
            min_radius=min_radius,
            min_pairs=min_pairs,
            min_variable_size=min_variable_size,
            max_variable_size=max_variable_size,
            min_constant_size=min_constant_size,
            substructure_pat=substructure_pat,
            pool=pool,
            explain=explain,
            )
    except analysis_algorithms.EvalError as err:
        sys.stderr.write("ERROR: %s\nExiting.\n" % (err,))
        raise SystemExit(1)

    transform_time = time.time()
    
    with outfile:
        result.write_products(
            outfile,
            field_names = (
#                "rule_environment_statistics_id",),
                "from_smiles", "to_smiles", "radius", "fingerprint", "rule_environment_id",
                "count", "avg", "std", "kurtosis", "skewness",
                "min", "q1", "median", "q3", "max", "paired_t", "p_value"),
            #column_aliases = {"from_smiles": "FROM"}, # use this to change the column name for a field
            include_empty=include_empty)

    output_time = time.time()
    
    if args.times:
        sys.stderr.write("Elapsed time (in seconds):\n")
        format_dt = get_time_delta_formatter(output_time - start_time)
        sys.stderr.write("  open database: %s\n" % format_dt(open_time - start_time))
        sys.stderr.write("  prepare query: %s\n" % format_dt(query_prep_time - open_time))
        sys.stderr.write("      transform: %s\n" % format_dt(transform_time - query_prep_time))
        sys.stderr.write("   write output: %s\n" % format_dt(output_time - transform_time))
        sys.stderr.write("         TOTAL = %s\n" % format_dt(output_time - start_time))
        
########################

