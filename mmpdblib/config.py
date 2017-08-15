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

from __future__ import absolute_import

#### Handle command-line arguments ####

from collections import OrderedDict
import argparse

from ._compat import basestring
from . import smarts_aliases

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
            raise argparse.ArgumentTypeError("could not parse %r as an integer: %s"
                                             % (term, err))
        
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

def parse_num_cuts_value(value):
    if value not in ("1", "2", "3"):
        raise argparse.ArgumentTypeError("must be '1', '2', or '3'")
    return int(value)
        

def parse_method_value(value):
    if value not in ("chiral",):
        raise argparse.ArgumentTypeError("must be 'chiral'")
    return value

class FragmentOptions(object):
    def __init__(self, max_heavies, max_rotatable_bonds,
                 rotatable_smarts, cut_smarts, num_cuts,
                 method, salt_remover):
        assert isinstance(max_heavies, int) or max_heavies is None, max_heavies
        self.max_heavies = max_heavies

        assert isinstance(max_rotatable_bonds, int) or max_rotatable_bonds is None, max_rotatable_bonds
        self.max_rotatable_bonds = max_rotatable_bonds

        assert isinstance(rotatable_smarts, str), rotatable_smarts
        self.rotatable_smarts = rotatable_smarts
        
        assert isinstance(cut_smarts, str), cut_smarts
        self.cut_smarts = cut_smarts

        assert num_cuts in (1, 2, 3), num_cuts
        self.num_cuts = num_cuts

        assert method in ("chiral",)
        self.method = method

        assert isinstance(salt_remover, basestring), salt_remover
        self.salt_remover = salt_remover

    def to_dict(self):
        d = OrderedDict()
        for name in ("max_heavies", "max_rotatable_bonds", "rotatable_smarts",
                     "cut_smarts", "num_cuts", "method", "salt_remover"):
            d[name] = getattr(self, name)
        return d
    
    def to_text_settings(self):
        def _none(x):
            return "none" if x is None else str(x)
        return (
            ("max_heavies", _none(self.max_heavies)),
            ("max_rotatable_bonds", _none(self.max_rotatable_bonds)),
            ("rotatable_smarts", self.rotatable_smarts),
            ("cut_smarts", self.cut_smarts),
            ("num_cuts", str(self.num_cuts)),
            ("method", self.method),
            ("salt_remover", self.salt_remover),
            )

DEFAULT_FRAGMENT_OPTIONS = FragmentOptions(
    max_heavies = 100,
    max_rotatable_bonds = 10,
    rotatable_smarts = "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]",
    cut_smarts = smarts_aliases.cut_smarts_aliases_by_name["default"].smarts,
    num_cuts = 3,
    method = "chiral",
    salt_remover = "<default>",
    )

def add_fragment_arguments(parser):
    p = parser
    OPTS = DEFAULT_FRAGMENT_OPTIONS
    p.add_argument("--max-heavies", type=parse_max_heavies_value,
                   metavar="N", default=None,
                   help="Maximum number of non-hydrogen atoms, or 'none' (default: %r)"
                         % (OPTS.max_heavies,))
    
    p.add_argument("--max-rotatable-bonds", type=parse_max_rotatable_bonds_value,
                   metavar="N", default=None,
                   help="Maximum number of rotatable bonds (default: %d)"
                         % (OPTS.max_rotatable_bonds,))
    
    p.add_argument("--rotatable-smarts", metavar="SMARTS",
                   help="SMARTS pattern to detect rotatable bonds (default: %r)"
                         % (OPTS.rotatable_smarts,))
    
    p.add_argument("--salt-remover", metavar="FILENAME",
                   help="File containing RDKit SaltRemover definitions. The default ('<default>') "
                        "uses RDKit's standard salt remover. Use '<none>' to not remove salts.")

    alias_names = ", ".join(repr(alias.name) for alias in smarts_aliases.cut_smarts_aliases)
    p.add_argument("--cut-smarts", metavar="SMARTS",
                   help="alternate SMARTS pattern to use for cutting (default: %r), or use one of: %s"
                   % (OPTS.cut_smarts, alias_names))
    
    p.add_argument("--num-cuts", choices=(1, 2, 3), type=parse_num_cuts_value,
                   help="number of cuts to use (default: %d)"
                   % (OPTS.num_cuts,))
    
    ## p.add_argument("--method", choices=("dalke", "hussain"), type=fragment_io.parse_method_value,
    ##                help="fragment canonicalization method to use (default: %s)"
    ##                % (FRAGMENT_OPTIONS.method,))

###### Index

parse_min_variable_heavies_value = nonnegative_int
parse_max_variable_heavies_value = positive_int_or_none

parse_max_variable_ratio_value = nonnegative_float
parse_min_variable_ratio_value = positive_float
parse_max_heavies_transf = nonnegative_int
parse_max_frac_trans = nonnegative_float

class IndexingOptions(object):
    def __init__(self,
                 min_variable_heavies, max_variable_heavies,
                 min_variable_ratio, max_variable_ratio,
                 max_heavies_transf, max_frac_trans
                 ):
        assert min_variable_heavies is None or min_variable_heavies >= 0, min_variable_heavies
        self.min_variable_heavies = min_variable_heavies

        assert (    (max_variable_heavies is None)
                 or (min_variable_heavies is None and max_variable_heavies >= 0)
                 or (min_variable_heavies <= max_variable_heavies)), max_variable_heavies
        self.max_variable_heavies = max_variable_heavies

        assert min_variable_ratio is None or 0.0 <= min_variable_ratio <= 1.0, min_variable_ratio
        self.min_variable_ratio = min_variable_ratio

        assert (   (max_variable_ratio is None)
                or (min_variable_ratio is not None and max_variable_ratio <= 1.0)
                or (min_variable_ratio <= max_variable_ratio <= 1.0) )
        self.max_variable_ratio = max_variable_ratio

        assert max_heavies_transf is None or max_heavies_transf >= 0, max_heavies_transf
        self.max_heavies_transf = max_heavies_transf

        assert max_frac_trans is None or max_frac_trans >= 0, max_heavies_transf
        self.max_frac_trans = max_frac_trans

        
DEFAULT_INDEX_OPTIONS = IndexingOptions(
    min_variable_heavies = None, # XXX can this be 0?
    max_variable_heavies = 10,
    min_variable_ratio = None,   # XXX can this be 0.0?
    max_variable_ratio = None,
    max_heavies_transf = None,
    max_frac_trans = None,  # XXX can this be 1.0?
    )

def add_index_options(parser):
    p = parser
    OPTS = DEFAULT_INDEX_OPTIONS
    p.add_argument("--min-variable-heavies", type=parse_min_variable_heavies_value,
                   metavar="N", default=OPTS.min_variable_heavies,
                    help="Minimum number of non-hydrogen atoms in the variable fragment.")
    p.add_argument("--max-variable-heavies", type=parse_max_variable_heavies_value,
                   default=DEFAULT_INDEX_OPTIONS.max_variable_heavies,
                   metavar="N",
                   help="Maximum number of non-hydrogen atoms in the variable fragment "
                        "(default: 10; for no maximum use 'none')")
    p.add_argument("--min-variable-ratio", type=parse_min_variable_ratio_value, default=None,
                   metavar="FLT",
                   help="Minimum ratio of variable fragment heavies to heavies in the (cleaned) structure")
    p.add_argument("--max-variable-ratio", type=parse_max_variable_ratio_value, default=None,
                   metavar="FLT",
                   help="Maximum ratio of variable fragment heavies to heavies in the (cleaned) structure")
    p.add_argument("--max-heavies-transf", type=parse_max_heavies_transf, default=None,
                   metavar="N",
                   help="Maximum difference in the number of heavies transfered in a transformation")
    p.add_argument("--max-frac-trans", type=parse_max_frac_trans, default=None,
                   metavar="FLT",
                   help="Maximum fraction of atoms taking part in a transformation")


class DEFAULT_RULE_SELECTION_OPTIONS:
    where = None
    score = None
    cutoff_list = (10, 5, 0)

def add_rule_selection_arguments(parser):
    OPTS = DEFAULT_RULE_SELECTION_OPTIONS
    parser.add_argument("--where", metavar="EXPR", default=OPTS.where,
                        help="select only rules for which the expression is true")
    parser.add_argument("--score", metavar="EXPR", default=OPTS.score,
                        help="use to break ties when multiple rules produce the same SMILES")
    parser.add_argument("--rule-selection-cutoffs", metavar="LIST", type=cutoff_list, default=OPTS.cutoff_list,
                        dest="cutoffs",
                        help="evaluate rule environments with the given minimum pair count. If multiple "
                             "counts are given, consider them in turn until there is a selected environment. "
                             "(default: '10,5,0')")
    
