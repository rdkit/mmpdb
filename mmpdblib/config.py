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

#### Handle command-line arguments ####

import argparse

from . import smarts_aliases
from . import fragment_types
from . import index_types

# Things to pass as the ArgumentParser argument's 'type'


def positive_int(value):
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def positive_int_or_none(value):
    if value == "none":
        return "none"
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'")
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'")
    return value


def positive_float(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive float")
    if not (value > 0.0):
        raise argparse.ArgumentTypeError("must be a positive float")
    return value


def nonnegative_float(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive float or zero")
    if not (value >= 0.0):
        raise argparse.ArgumentTypeError("must be a positive float or zero")
    return value


def nonnegative_int(value):
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer or zero")
    if not (value >= 0):
        raise argparse.ArgumentTypeError("must be a positive integer or zero")
    return value


def cutoff_list(value_s):
    prev = None
    values = []
    for term in value_s.split(","):
        try:
            value = int(term)
        except ValueError as err:
            raise argparse.ArgumentTypeError("could not parse %r as an integer: %s" % (term, err))

        if value < 0:
            raise argparse.ArgumentTypeError("threshold values must be non-negative")

        if prev is not None and prev <= value:
            raise argparse.ArgumentTypeError("threshold values must be in decreasing order")
        prev = value

        values.append(value)

    if not values:  # Let people specify ""
        return [0]

    return values


#### Fragment

parse_max_heavies_value = positive_int_or_none
parse_max_rotatable_bonds_value = positive_int_or_none
parse_min_heavies_per_const_frag_value = nonnegative_int
parse_min_heavies_total_const_frag_value = nonnegative_int


def parse_num_cuts_value(value):
    if value not in ("1", "2", "3"):
        raise argparse.ArgumentTypeError("must be '1', '2', or '3'")
    return int(value)


def parse_method_value(value):
    if value not in ("chiral",):
        raise argparse.ArgumentTypeError("must be 'chiral'")
    return value


DEFAULT_FRAGMENT_OPTIONS = fragment_types.FragmentOptions(
    max_heavies=100,
    max_rotatable_bonds=10,
    rotatable_smarts="[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]",
    cut_smarts=smarts_aliases.cut_smarts_aliases_by_name["default"].smarts,
    num_cuts=3,
    method="chiral",
    salt_remover="<default>",
    min_heavies_per_const_frag=0,
    min_heavies_total_const_frag=0,
    max_up_enumerations=1000,
)


###### Index

parse_min_variable_heavies_value = nonnegative_int
parse_max_variable_heavies_value = positive_int_or_none

parse_max_variable_ratio_value = nonnegative_float
parse_min_variable_ratio_value = positive_float
parse_max_heavies_transf = nonnegative_int
parse_max_frac_trans = nonnegative_float
parse_max_radius = nonnegative_int


DEFAULT_INDEX_OPTIONS = index_types.IndexOptions(
    min_variable_heavies=None,  # XXX can this be 0?
    max_variable_heavies=10,
    min_variable_ratio=None,  # XXX can this be 0.0?
    max_variable_ratio=None,
    max_heavies_transf=None,
    max_frac_trans=None,  # XXX can this be 1.0?,
    min_radius=0,
    max_radius=5,
    symmetric=False,
    smallest_transformation_only=False,
)


def add_index_options(parser):
    p = parser
    OPTS = DEFAULT_INDEX_OPTIONS
    p.add_argument(
        "--min-variable-heavies",
        type=parse_min_variable_heavies_value,
        metavar="N",
        default=OPTS.min_variable_heavies,
        help="Minimum number of non-hydrogen atoms in the variable fragment.",
    )
    p.add_argument(
        "--max-variable-heavies",
        type=parse_max_variable_heavies_value,
        default=DEFAULT_INDEX_OPTIONS.max_variable_heavies,
        metavar="N",
        help="Maximum number of non-hydrogen atoms in the variable fragment "
        "(default: 10; for no maximum use 'none')",
    )
    p.add_argument(
        "--min-variable-ratio",
        type=parse_min_variable_ratio_value,
        default=None,
        metavar="FLT",
        help="Minimum ratio of variable fragment heavies to heavies in the (cleaned) structure",
    )
    p.add_argument(
        "--max-variable-ratio",
        type=parse_max_variable_ratio_value,
        default=None,
        metavar="FLT",
        help="Maximum ratio of variable fragment heavies to heavies in the (cleaned) structure",
    )
    p.add_argument(
        "--max-heavies-transf",
        type=parse_max_heavies_transf,
        default=None,
        metavar="N",
        help="Maximum difference in the number of heavies transfered in a transformation",
    )
    p.add_argument(
        "--max-frac-trans",
        type=parse_max_frac_trans,
        default=None,
        metavar="FLT",
        help="Maximum fraction of atoms taking part in a transformation",
    )
    p.add_argument(
        "--max-radius",
        type=parse_max_radius,
        default=DEFAULT_INDEX_OPTIONS.max_radius,
        metavar="N",
        help="Maximum Environment Radius to be indexed in the MMPDB database",
    )


class RuleSelectionOptions:
    def __init__(self, where, score, cutoff_list):
        self.where = where
        self.score = score
        self.cutoff_list = cutoff_list

    def get_rule_selection_function(self):
        from . import analysis_algorithms

        # generate the sort key for the different cutoffs
        score_function = self.score
        if score_function is None:
            score_function = analysis_algorithms.default_score_function

        rule_key_function = analysis_algorithms.ComputeRuleKey(
            score_function=score_function,
            cutoffs=self.cutoff_list,
        )
        # wrap it together into a selection function
        return analysis_algorithms.RuleSelectionFunction(
            self.where,
            rule_key_function,
        )


DEFAULT_RULE_SELECTION_OPTIONS = RuleSelectionOptions(
    where=None,
    score=None,
    cutoff_list=(10, 5, 0),
)
