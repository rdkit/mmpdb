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

from __future__ import print_function, absolute_import, division

import sys
import os

from . import command_support
from . import index_algorithm
from . import fragment_io
from . import fragment_types
from . import properties_io
from ._compat import open_universal

try:
    import psutil
except ImportError:
    psutil = None

def get_fragment_filter_from_args(parser, args):
    min_variable_heavies = args.min_variable_heavies
    max_variable_heavies = args.max_variable_heavies
    min_variable_ratio = args.min_variable_ratio
    max_variable_ratio = args.max_variable_ratio

    if max_variable_heavies == "none": # Push this into the argument parser?
        max_variable_heavies = None

    if min_variable_ratio is not None and max_variable_ratio is not None:
        if min_variable_ratio > max_variable_ratio:
            parser.error("--min-variable-ratio must not be larger than --max-variable-ratio")

    if min_variable_heavies is not None and max_variable_heavies is not None:
        if min_variable_heavies > max_variable_heavies:
            parser.error("--min-variable-heavies must not be larger than --max-variable-heavies")
        
        
    filters = []
    if min_variable_heavies is not None:
        filters.append(index_algorithm.MinVariableHeaviesFilter(min_variable_heavies))
    if max_variable_heavies is not None:
        filters.append(index_algorithm.MaxVariableHeaviesFilter(max_variable_heavies))
    if min_variable_ratio is not None:
        filters.append(index_algorithm.MinVariableRatioFilter(min_variable_ratio))
    if max_variable_ratio is not None:
        filters.append(index_algorithm.MaxVariableRatioFilter(max_variable_ratio))

    if not filters:
        # It's easier to have 0 filters than to make a special do-nothing filter.
        return index_algorithm.MultipleFilters([])
    elif len(filters) == 1:
        return filters[0]
    else:
        return index_algorithm.MultipleFilters(filters)

if psutil is None:
    def get_memory_use():
        return 0
else:
    _process = psutil.Process(os.getpid())
    def get_memory_use():
        info = _process.memory_info()
        return info.rss #  or info.vms?

def human_memory(n):
    if n < 1024:
        return "%d B" % (n,)
    for unit, denom in (("KB", 1024), ("MB", 1024**2),
                         ("GB", 1024**3), ("PB", 1024**4)):
        f = n/denom
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
    
def index_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    if args.title:
        title = args.title
    else:
        title="MMPs from %r" % (args.fragment_filename,)

    fragment_filter = get_fragment_filter_from_args(parser, args)

    report_memory = False
    if args.memory:
        if psutil is None:
            sys.stderr.write("WARNING: Cannot report memory because the 'psutil' module is not available.\n")
        else:
            report_memory = True

    start_properties_memory = get_memory_use()
    if args.properties is None:
        selected_ids = properties = None
    else:
        try:
            properties_file = open_universal(args.properties)
        except IOError as err:
            parser.error("Cannot open --properties file: %s" % (err,))

        try:
            with properties_file:
                properties = properties_io.load_properties(properties_file, reporter)
        except ValueError as err:
            parser.error("Problem reading --properties file %r: %s"
                         % (args.properties, err))
            
        selected_ids = set(properties.get_ids())
    end_properties_memory = get_memory_use()
    

    if ((args.out is None or args.out == "mmpdb") and args.output is None):
        # Use the filename based on the fragments filename
        fragment_filename = args.fragment_filename
        if fragment_filename is None:
            parser.error("The '--out mmpdb' format requires a filename when reading from stdin.")

        # replace the extension (if any) with ".mmpdb"
        args.output = os.path.splitext(fragment_filename)[0] + ".mmpdb"
        reporter.warning("No --output filename specified. Saving to %r." % (args.output,))

    #reporter.report("Using fragment filters: %s" % (fragment_filter.get_args(),))

    fragment_io.suggest_faster_json(reporter)

    start_fragment_index_memory = get_memory_use()
    try:
        fragment_reader = fragment_io.read_fragment_records(args.fragment_filename)
    except fragment_types.FragmentFormatError as err:
        parser.error(str(err))
        
    with fragment_reader:
        with reporter.progress(fragment_reader, "Loaded fragment record") as report_fragment_reader:
            fragment_index = index_algorithm.load_fragment_index(report_fragment_reader, fragment_filter, selected_ids)

    index_options = fragment_filter.get_options()
    index_options.update(
        symmetric = args.symmetric,
        max_heavies_transf = args.max_heavies_transf,
        max_frac_trans = args.max_frac_trans,
        )
    index_options = index_algorithm.IndexOptions(**index_options)
                           
    start_mmp_memory = get_memory_use()
    environment_cache = index_algorithm.EnvironmentCache()
    pairs = index_algorithm.find_matched_molecular_pairs(
        fragment_index, index_options,
        reporter=reporter)


    with index_algorithm.open_mmpa_writer(args.output, format=args.out,
                                   title=title,
                                   fragment_options=fragment_reader.options,
                                   fragment_index=fragment_index,
                                   index_options=index_options,
                                   properties=properties,
                                   environment_cache=environment_cache,
                                   ) as pair_writer:
        pair_writer.start()
        pair_writer.write_matched_molecule_pairs(pairs)
        end_mmp_memory = get_memory_use()
        pair_writer.end(reporter)
        end_memory = get_memory_use()
        
    if report_memory:
        #print(start_fragment_index_memory, start_mmp_memory, end_mmp_memory, end_memory)
        sys.stderr.write("#pairs: %d Memory (RSS) total: %s properties: %s fragments: %s indexing: %s\n" % (
                         pair_writer.num_pairs, 
                         human_memory(max(end_mmp_memory, end_memory)),
                         human_memory(end_properties_memory-start_properties_memory),
                         human_memory(start_mmp_memory - start_fragment_index_memory),
                         human_memory(end_mmp_memory - start_mmp_memory)))
        

