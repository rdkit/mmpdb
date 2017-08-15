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

from . import reporters
from . import fileio

def get_reporter(is_quiet):
    if is_quiet:
        return reporters.get_reporter("quiet")
    else:
        return reporters.get_reporter("verbose")

def explain(msg, *args):
    full_msg = (msg % args) + "\n"
    sys.stderr.write(full_msg)

def no_explain(msg, *args):
    pass

def get_explain(use_explain, reporter=None):
    if use_explain:
        if reporter is None:
            return explain
        else:
            return reporter.explain
    return no_explain


def get_property_names_or_error(parser, args, dataset):
    property_names = []
    seen = set()

    if args.property:
        if args.no_properties:
            parser.error("Cannot specify --property and --no-properties")
        if args.all_properties:
            parser.error("Cannot specify --property and --all-properties")
    
    known_names = dataset.get_property_names()
    if not args.property:
        if args.no_properties:
            return []
        # Use all of the properties
        return known_names

    known_names = set(known_names)
    for name in args.property:
        # If a property is specified multiple times, only use the first one.
        if name in seen:
            continue
        seen.add(name)

        if name not in known_names:
            parser.error("--property %r is not present in the database"
                         % (name,))
        property_names.append(name)
    return property_names
        


def get_property_name_or_error(parser, args, db, dataset):
    assert args.property is not None # the --property must be 'required'
    property_name = dataset.get_property_name_with_name(args.property)
    if property_name is not None:
        return property_name
    
    # Not found. Provide a little extra in the error
    names = [property_name.name for property_name in dataset.get_property_names_table().values()]
    if names:
        extra = "(Available properties: %s)" % (repr(names)[1:-1],)
    else:
        extra = "(No properties available)"
    parser.error(
        "--property %r not found in this dataset\n%s"
        % (args.property, extra))

