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

from __future__ import print_function, absolute_import

from collections import defaultdict
import re

from rdkit import Chem

from . import command_support
from . import do_fragment
from . import fragment_algorithm
from . import index_algorithm
from . import environment
from . import smiles_syntax
from . import schema
from .config import DEFAULT_RULE_SELECTION_OPTIONS
from . import _compat

####

class EvalError(Exception):
    pass

## Some helper object which are always ranked at the top or the bottom
class SmallestScore(object):
    def __lt__(self, other):
        if isinstance(other, SmallestScore):
            return False
        return True
    def __gt__(self, other):
        return False
    def __eq__(self, other):
        return isinstance(other, SmallestScore)
_smallest_score = SmallestScore()

class NInfinity(object):
    def __lt__(self, other):
        if isinstance(other, NInfinity):
            return False
        return True
    def __gt__(self, other):
        return False
    def __eq__(self, other):
        return isinstance(other, NInfinity)
    def __neg__(self):
        return infinity
    def __pos__(self):
        return self
ninfinity = NInfinity()  # represent negative infinity score

class Infinity(object):
    def __lt__(self, other):
        return False
    def __gt__(self, other):
        if isinstance(other, Infinity):
            return False
        return True
    def __eq__(self, other):
        return isinstance(other, Infinity)
    def __neg__(self):
        return ninfinity
    def __pos__(self):
        return self
infinity = Infinity()  # represent positive infinity score

##

# Check that the eval'ed code only references known variable names.
# (The goal is to provide some helpful errors in case of mistakes.
# Python's eval is *not* secure. There are ways to reach global
# objects without going through a variable name lookup.)
def check_eval_names(code_obj, extra_names):
    known_names = set(schema.PropertyRule.__slots__)
    known_names.update(extra_names)

    for name in code_obj.co_names:
        if name not in known_names:
            raise ValueError("unsupported variable name %r"
                             % (name,))
    

def get_where_function(where_expr=None):
    if where_expr is None:
        return None

    # Check that it's an eval-able term
    where_expr = where_expr.strip()
    try:
        where_code = compile(where_expr, "--where", "eval")
    except Exception as err:
        raise ValueError("Cannot parse: %s" % (err,))
    
    check_eval_names(where_code, ["__builtins__", "None", "True", "False"])

    def where_function(property_rule):
        d = property_rule.to_dict()
        d.update({
            "__builtins__": None,
            "None": None,
            "True": True,
            "False": False,
            })
        
        try:
            result = eval(where_code, d)
        except Exception as err:
            raise EvalError("Could not evaluate 'where' expression %r: %s" % (where_expr, err))

        return result
    
    return where_function

def default_score_function(property_rule):
    
    # Break ties on the minimum standard deviation
    std = property_rule.std
    if std is None: # only one element, so effectively infinite stddev
        std = ninfinity
    else:
        std = -std  # smallest std is best, so invert

    # Break ties on the rule environment with the highest radius
    radius = property_rule.radius

    # Break ties on the rule environment with the highest
    # count of heavy atoms in the LHS fragment
    lhs_num_heavies = property_rule.from_num_heavies

    # Break ties arbitrarily on the LHS SMILES
    lhs_smiles = property_rule.from_smiles

    return (std, radius, lhs_num_heavies, lhs_smiles)


def get_score_function(score_expr=None):
    if score_expr is None:
        return default_score_function

    # .strip() because the command-line "--score '-std'" is interpreted as --score
    # with a missing value, and the option -s/-t/-d.
    # The workaround is: "--score ' -std'"
    score_expr = score_expr.strip()
    try:
        score_code = compile(score_expr, "--score", "eval")
    except Exception as err:
        raise ValueError("Cannot parse: %s" % (err,))

    check_eval_names(score_code,
                     ["__builtins__", "None", "True", "False", "inf", "ninf"])

    def score_function(property_rule):
        d = property_rule.to_dict()
        d.update({
            "__builtins__": None,
            "None": None,
            "inf": infinity,
            "ninf": ninfinity,
            "True": True,
            "False": False,
            })
        try:
            result = eval(score_code, d)
        except Exception as err:
            raise EvalError("Could not evaluate --score %r: %s" % (score_expr, err))

        return result

    return score_function

class ComputeRuleKey(object):
    def __init__(self, score_function, cutoffs):
        self.score_function = score_function
        self.cutoffs = cutoffs

    def __call__(self, property_rule):
        score = self.score_function(property_rule)
        key = []
        
        for cutoff in self.cutoffs:
            if property_rule.count < cutoff:
                # Force this to be the smallest value and try the next cutoff
                key.append( (0, _smallest_score) )
            else:
                key.append( (1, score) )

        return tuple(key)


class RuleSelectionFunction(object):
    def __init__(self, where_function, rule_key_function):
        self.where_function = where_function
        self.rule_key_function = rule_key_function

    def __call__(self, property_rules, explain=command_support.no_explain):
        if not property_rules:
            explain("    No rule environment statistics has enough pairs.")
            return None

        # Order the list so the --explain is easier to understand.
        def by_rule_id_then_radius(property_rule):
            return (property_rule.rule_id, property_rule.radius)
        property_rules = sorted(property_rules, key=by_rule_id_then_radius)
        
        # Apply the --where clause, if specified.
        if self.where_function is not None:
            property_rules = self.apply_where_function(property_rules, explain)
            if property_rules is None:
                return None

        # Find the highest-ranked rule, if any of the rules meet the minimum cutoff.
        result = self.select_max_rule(property_rules, explain)
        if result is None:
            return None
        property_rule, property_rules = result

        # Report what happened.
        _explain_choice(property_rule, property_rules, explain)
        
        return property_rule
        

    def apply_where_function(self, property_rules, explain):
        new_list = []
        for property_rule in property_rules:
            # Always do the test, even if there is only one term.
            if self.where_function(property_rule):
                new_list.append(property_rule)
                explain("    Rule %s (#%d) passed the --where test",
                        property_rule.smirks, property_rule.rule_id)
            else:
                explain("    Rule %s (#%d) did not pass the --where test",
                        property_rule.smirks, property_rule.rule_id)

        if not new_list:
            explain("    No rule environment statistics passed the --where test.")
            return None
            
        return new_list

    def select_max_rule(self, property_rules, explain):
        # Select the best option
        property_rules = sorted(property_rules,
                                          key=self.rule_key_function, reverse=True)
        property_rule = property_rules[0]

        # Make sure it's actually usable. Otherwise it could be that
        # none of the rules have enough pairs for the cutoffs.
        key = self.rule_key_function(property_rule)
        for in_use, score in key:
            if in_use:
                return property_rule, property_rules

        # No rules meet any of the cutoffs. Ignore.
        explain("    No rule environment statistics meets the minimum cutoff")
        return None

def _explain_choice(property_rule, property_rules, explain):
    if len(property_rules) == 1:
        explain("Selected rule %s (rule %d).",
                property_rule.smirks, property_rule.rule_id)
    else:
        explain("Selected rule %s (#%d). The %d possibilities were:",
                property_rule.smirks, property_rule.rule_id, len(property_rules))
        for property_rule in property_rules:
            explain("   %s (rule %d) radius: %d count: %d std: %s",
                    property_rule.smirks, property_rule.rule_id, property_rule.radius,
                    property_rule.count, str(property_rule.std) if property_rule.std is not None else "n/a",
                    )
            ## for q in self.rule_key_function(stats):
            ##     explain("     %r", q)
                        
    
def get_rule_selection_function(where_expr, score_expr, cutoffs):
    # '--where': select which rule environments to use
    where_function = get_where_function(where_expr)
    # '--score': rank the rule environments
    score_function = get_score_function(score_expr)
    # generate the sort key for the different cutoffs
    rule_key_function = ComputeRuleKey(score_function=score_function, cutoffs=cutoffs)
    # wrap it together into a selection function
    return RuleSelectionFunction(where_function, rule_key_function)

def get_rule_selection_function_from_args(parser, args):
    try:
        where_function = get_where_function(args.where)
    except ValueError as err:
        parser.error("%s, in --where %r" % (err, args.where))
        
    try:
        score_function = get_score_function(args.score)
    except ValueError as err:
        parser.error("%s, in --score %r" % (err, args.score))
        
    rule_key_function = ComputeRuleKey(score_function=score_function, cutoffs=args.cutoffs)
    
    return RuleSelectionFunction(where_function, rule_key_function)
    

default_rule_selection_function = get_rule_selection_function(
    where_expr=DEFAULT_RULE_SELECTION_OPTIONS.where,
    score_expr=DEFAULT_RULE_SELECTION_OPTIONS.score,
    cutoffs=DEFAULT_RULE_SELECTION_OPTIONS.cutoff_list)



        



###### Predict

def _get_tool(klass, dataset, rule_selection_function):
    cursor = dataset.get_cursor()
    property_name_to_id = dataset.get_property_names_table(cursor)
    fragment_options = dataset.get_fragment_options(cursor)
    fragment_filter = do_fragment.get_fragment_filter(fragment_options)
    return klass(
        dataset = dataset,
        property_name_to_id = property_name_to_id,
        fragment_options = fragment_options,
        fragment_filter = fragment_filter,
        rule_selection_function = rule_selection_function,
        )

class Tool(object):
    def __init__(self,
                 dataset, property_name_to_id,
                 fragment_options, fragment_filter, rule_selection_function):
        self.dataset = dataset
        self.property_name_to_id = property_name_to_id
        self.fragment_options = fragment_options
        self.fragment_filter = fragment_filter
        self.rule_selection_function = rule_selection_function

    def get_property_names(self):
        return list(self.property_name_to_id)
    
    def is_available_property_name(self, property_name):
        return property_name in self.property_name_to_id

    
def get_predict_tool(
        dataset,
        rule_selection_function = default_rule_selection_function
        ):
    return _get_tool(PredictTool, dataset, rule_selection_function)


class PredictTool(Tool):
    def fragment_predict_smiles(self, smiles):
        if "[H]" in smiles:
            smiles_record = do_fragment.make_hydrogen_fragment_record("query", smiles, self.fragment_filter)
            return smiles_record

        smiles_record = do_fragment.make_fragment_record_from_smiles(smiles, self.fragment_filter)
        if smiles_record.errmsg is not None:
            return smiles_record

        # Include hydrogen fragments
        hydrogen_fragments = fragment_algorithm.get_hydrogen_fragmentations(
            smiles_record.normalized_smiles, smiles_record.num_normalized_heavies)

        smiles_record.fragments.extend(hydrogen_fragments)
        return smiles_record

    def fragment_reference_smiles(self, smiles):
        smiles_record = do_fragment.make_fragment_record_from_smiles(smiles, self.fragment_filter)
        if smiles_record.errmsg is not None:
            return smiles_record

        # Include hydrogen fragments
        hydrogen_fragments = fragment_algorithm.get_hydrogen_fragmentations(
            smiles_record.normalized_smiles, smiles_record.num_normalized_heavies)

        smiles_record.fragments.extend(hydrogen_fragments)
        return smiles_record

    def predict(self, reference_fragments, smiles_fragments,
                property_name, explain=None, relabel_cache=None):
        property_name_id = self.property_name_to_id.get(property_name, None)
        if property_name_id is None:
                raise ValueError("Property %r not found")
        cursor = self.dataset.get_cursor()
        return make_prediction(
            self.dataset, reference_fragments, smiles_fragments,
            property_name, property_name_id,
            rule_selection_function=self.rule_selection_function,
            relabel_cache=relabel_cache, cursor=cursor, explain=explain)
    
#####

# Get the attachment numbers from a SMILES string.
# For example: "C[*:1]N([*:3])P[*:2]" returns "132".
# XXX or should this return inverse?
_attachment_pattern = re.compile(":([123])")  
def get_attachment_order(smiles):
    order = "".join(_attachment_pattern.findall(smiles))
    return order

    
def make_prediction(
        dataset, reference_fragments, smiles_fragments,
        property_name, property_name_id,
        rule_selection_function,
        explain=None,
        relabel_cache=None,
        cursor=None,
        ):
    if explain is None:
        explain = command_support.no_explain
        using_explain = False
    else:
        using_explain = True
    
    # Used to help speed up the canonical SMIRKS.
    # It's passed in as an argument because it might be shared across multiple
    # 'predict' and 'transform' calls/datasets, e.g., for a hypothetical web
    # service. Note: that system would need to implement occasional pruning
    # because the current cache is unbounded.
    relabel_cache = index_algorithm.RelabelCache()

    if cursor is None:
        cursor = dataset.get_cursor()

    # Find the variable fragments for each constant fragment in --reference
    constant_to_reference_fragments = defaultdict(list)
    for frag in reference_fragments:
        explain("--reference fragments: constant %s variable %s", frag.constant_smiles, frag.variable_smiles)
        constant_to_reference_fragments[frag.constant_smiles].append(frag)
        
    # Find the variable fragments for each constant fragment in --smiles
    constant_to_smiles_fragments = defaultdict(list)
    for frag in smiles_fragments:
        explain("--smiles fragments: constant %s variable %s", frag.constant_smiles, frag.variable_smiles)
        constant_to_smiles_fragments[frag.constant_smiles].append(frag)

    # Which constants are common to both --reference and --smiles?
    common_constants = set(constant_to_reference_fragments) & set(constant_to_smiles_fragments)

        
    # Figure out which rules to consider for each constant
    unique_property_rules = {}
    for constant_smiles in common_constants:
        explain("Found match for constant %s (#variable fragments in the --reference and --smiles: %dx%d)",
                constant_smiles,
                len(constant_to_reference_fragments[constant_smiles]),
                len(constant_to_smiles_fragments[constant_smiles]))

        # Get the environment fingerprints around each attachment point
        all_center_fps = environment.compute_constant_center_fingerprints(constant_smiles)

        # For every possible rule from LHS to RHS ...
        for from_frag in constant_to_reference_fragments[constant_smiles]:
            for to_frag in constant_to_smiles_fragments[constant_smiles]:
                assert from_frag.constant_symmetry_class == to_frag.constant_symmetry_class, (
                    constant_smiles, from_frag.variable_smiles, to_frag.variable_smiles,
                    from_frag.constant_symmetry_class, to_frag.constant_symmetry_class)
                constant_symmetry_class = to_frag.constant_symmetry_class

                # Skip the idempotent transformation. (Can occur if --smiles and --reference are identical.)
                if from_frag.variable_smiles == to_frag.variable_smiles:
                    continue

                # Figure out the canonical SMIRKS and canonical form of the constant.
                smirks, ordered_constant = index_algorithm.cansmirks(
                    from_frag.num_cuts,
                    from_frag.variable_smiles, from_frag.variable_symmetry_class, from_frag.attachment_order,
                    constant_smiles, from_frag.constant_symmetry_class,
                    to_frag.variable_smiles, to_frag.variable_symmetry_class, to_frag.attachment_order,
                    relabel_cache)
                explain("  SMIRKS %s with constant assignment %s", smirks, ordered_constant)

                # Use the ordered_constant to get the attachment order for the constant
                from_smiles, ignore, to_smiles = smirks.partition(">>")
                assert ignore == ">>", smirks
                new_constant_order = get_attachment_order(ordered_constant)

                # From that, compute the environment fingerprint(s)
                possible_envs = environment.get_all_possible_fingerprints(
                    all_center_fps, constant_symmetry_class, new_constant_order)
                #explain("  Environment fingerprints: %s", possible_envs)

                # Turn those into identifiers because it helps with the
                # iter_selected_property_rules() performance. I originally
                # tried to use the fingerprint string directly, but the SQL
                # query optimizer didn't optimize it very well.
                
                allowed_fingerprint_ids = dataset.get_fingerprint_ids(possible_envs, cursor)

                # Find the rules for the given LHS/RHS and property
                property_rules = dataset.iter_selected_property_rules(
                    from_smiles, to_smiles, property_name_id, cursor)
                
                # Select only those with the right environment fingerprints.
                property_rules = [rule for rule in property_rules
                                    if rule.fingerprint_id in allowed_fingerprint_ids]
                    
                if using_explain:
                    # Give a more in-depth description of the selected rules.
                    if not property_rules:
                        explain("  No rule found for %s>>%s or its inverse",
                                from_smiles, to_smiles)
                    else:
                        for rule in property_rules:
                            explain("  Rule %s (#%d) environment %d (radius %d): count %d avg %g",
                                    rule.smirks, rule.rule_id, rule.rule_environment_id,
                                    rule.radius, rule.count, rule.avg)
                            
                # De-duplicate the rules. This can happen if the same rule can be
                # used multiple times in a structure.
                for rule in property_rules:
                    unique_property_rules[ (rule.rule_environment_id, rule.is_reversed) ] = rule

    # Get all of the unique environments (includes all of the radii)
    property_rules = list(unique_property_rules.values())

    # Sort them so output is human-understandable
    def get_sort_key(rule):
        return (rule.rule_id, rule.radius, -rule.from_num_heavies)
    property_rules.sort(key = get_sort_key)

    # Figure out which ones are used in both directions
    directions = defaultdict(set)
    for rule in property_rules:
        directions[rule.rule_environment_id].add(rule.is_reversed)
    for rule in property_rules:
        rule.is_bidirectional = len(directions[rule.rule_environment_id]) > 1
    
    # Apply the rule selection function to get *the* property rule
    property_rule = rule_selection_function(property_rules, explain)

    return PredictResult(dataset, property_rule, property_rules,
                         property_name, property_name_id)


# Help generate table output
def _get_column_names(field_names, column_aliases):
    if not column_aliases:
        return field_names
    return [column_aliases.get(field_name, field_name) for field_name in field_names]

_global_column_formatters = {
    "rule_id": "%d",
    "rule_environment_id": "%d",
    "rule_environment_statistics_id": "%d",
    "radius": "%d",
    "from_smiles": "%s",
    "to_smiles": "%s",
    "fingerprint": "%s",
    "count": "%d",
    "avg": "%.5g",
    "std": "%.5g",
    "kurtosis": "%.5g",
    "skewness": "%.5g",
    "min": "%.5g",
    "q1": "%.5g",
    "median": "%.5g",
    "q3": "%.5g",
    "max": "%.5g",
    "paired_t": "%.5g",
    "p_value": "%.5g",

    "pair_id": "%d",
    "lhs_smiles": "%s",
    "rhs_smiles": "%s",
    "lhs_public_id": "%s",
    "rhs_public_id": "%s",
    "lhs_value": "%.5g",
    "rhs_value": "%.5g",
    "delta": "%.5g",
    }

def _format_object_attributes(obj, field_names, missing_value, formatters):
    value_strs = []
    for field_name in field_names:
        value = getattr(obj, field_name)
        if value is None:
            value_str = missing_value
        else:
            column_formatter = formatters[field_name]
            value_str = column_formatter % (value,)
        value_strs.append(value_str)
    return value_strs

def _format_object_attributes_to_string(obj, field_name, missing_value, formatters):
    return "\t".join(_format_object_attributes(
        obj, field_name, missing_value, formatters))

# Store the results of the transform prediction, and provide ways to get/save the data.
class PredictResult(object):
    missing_value = ""
    property_rule_field_names = (
        "rule_id", "rule_environment_id", "radius", "fingerprint", "from_smiles", "to_smiles",
        "rule_environment_statistics_id", "count", "avg", "std", "kurtosis", "skewness",
        "min", "q1", "median", "q3", "max", "paired_t", "p_value")

    property_rule_pair_field_names = (
        "rule_id", "rule_environment_id", "radius", "fingerprint", "pair_id",
        "lhs_smiles", "rhs_smiles", "lhs_public_id", "rhs_public_id",
        "lhs_value", "rhs_value", "delta")
        
    column_formatters = _global_column_formatters.copy()
    
    def __init__(self, dataset, property_rule, property_rules,
                     property_name, property_name_id):
        self.dataset = dataset
        self.property_rule = property_rule  # the selected rule
        self.property_rules = property_rules # all of the rules
        self.property_name = property_name
        self.property_name_id = property_name_id

    def format_property_rule(self, property_rule, field_names=None,
                             show_bidirectional_flag=True):
        if field_names is None:
            field_names = PredictResult.property_rule_field_names
    
        column_values = _format_object_attributes(
            property_rule, field_names, self.missing_value, self.column_formatters)
        
        if show_bidirectional_flag and property_rule.is_bidirectional:
            for i, field_name in field_names:
                if field_name == "avg":
                    column_values[i] += "*"
        return "\t".join(column_values)
        
    def write_property_rules(self, property_rules_file, field_names=None,
                             column_aliases=None,
                             show_bidirectional_flag=True):
        if field_names is None:
            field_names = PredictResult.property_rule_field_names
        column_names = _get_column_names(field_names, column_aliases)

        # write the header
        property_rules_file.write("\t".join(column_names) + "\n")

        # write a row for each rule
        for rule in self.property_rules:
            line = self.format_property_rule(
                rule, field_names, show_bidirectional_flag)
            property_rules_file.write(line + "\n")

    def get_property_rule_pairs(self, cursor=None):
        if cursor is None:
            cursor = self.dataset.get_cursor()
        for property_rule in self.property_rules:
            property_rule_pairs = self.dataset.get_property_rule_pairs(
                property_rule, self.property_name_id, cursor)
            for property_rule_pair in property_rule_pairs:
                yield property_rule_pair
            
    def format_property_rule_pair(self, rule_pair, field_names):
        if field_names is None:
            field_names = PredictResult.pair_field_names
        return _format_object_attributes_to_string(
            rule_pair, field_names, self.missing_value, self.column_formatters)
            
    def write_property_rule_pairs(self, rule_pairs_file, field_names=None,
                                  column_aliases=None):
        if field_names is None:
            field_names = PredictResult.property_rule_pair_field_names
        column_names = _get_column_names(field_names, column_aliases)

        # write the header
        rule_pairs_file.write("\t".join(column_names) + "\n")

        # write a row for each pair
        for rule_pair in self.get_property_rule_pairs():
            line = self.format_property_rule_pair(rule_pair, field_names)
            rule_pairs_file.write(line + "\n")
    

def test_predict():
    # --smiles 'c1ccccc1C(=O)N(C)C' --reference 'c1ccccc1C(=O)NC' --property logD mdo_mmp_logD.db --save-details
    import sys
    from . import dbutils
    db = dbutils.open_database("e.mmpdb")
    dataset = db.get_dataset()
    predict_tool = get_predict_tool(dataset)
    reference_record = predict_tool.fragment_reference_smiles("c1ccccc1C(=O)NC")
    if reference_record.errmsg is not None:
        raise ValueError("Cannot parse reference: %s" % (reference_record.errmsg,))
    
    predict_record = predict_tool.fragment_predict_smiles("c1ccccc1C(=O)N(C)C")
    if predict_record.errmsg is not None:
        raise ValueError("Cannot parse smiles: %s" % (predict_record.errmsg,))

    predict_results = predict_tool.predict(reference_record.fragments,
                                           predict_record.fragments, "MP")
    predict_results.write_property_rules(sys.stdout)
    predict_results.write_property_rule_pairs(sys.stdout)

###### Transform

def get_transform_tool(
        dataset,
        rule_selection_function = default_rule_selection_function
        ):
    return _get_tool(TransformTool, dataset, rule_selection_function)
    

class TransformTool(Tool):
    def fragment_transform_smiles(self, smiles):
        # Figure out how I'm going to fragment the input --smiles
        if "[H]" in smiles:
            # User-specified transform location
            record = do_fragment.make_hydrogen_fragment_record("query", smiles, self.fragment_filter)
        else:
            record = do_fragment.make_fragment_record_from_smiles(smiles, self.fragment_filter)
        return record
    
    def transform(self, transform_fragments, property_names,
                  min_radius=0, min_pairs=0, min_variable_size=0, min_constant_size=0,
                  substructure_pat=None,
                  pool=None,
                  explain=None):
        property_info_list = []
        for property_name in property_names:
            try:
                property_name_id = self.property_name_to_id[property_name]
            except KeyError:
                raise ValueError("Dataset does not contain property %r"
                                 % (property_name,))
            property_info_list.append( (property_name_id, property_name) )

        cursor = self.dataset.get_cursor()
        return make_transform(
                self.dataset, transform_fragments, property_info_list,
                rule_selection_function=self.rule_selection_function,
                substructure_pat=substructure_pat,
                min_radius=min_radius, min_pairs=min_pairs,
                min_variable_size=min_variable_size, min_constant_size=min_constant_size,
                pool=pool,
                cursor=cursor, explain=explain)


# Enumerate all of the ways that the canonical unlabeled SMILES
# might be turned into a non-canonical labeled SMILES.

_bracket_wildcard_pat = re.compile(re.escape("[*]"))
_organic_wildcard_pat = re.compile(re.escape("*"))

def enumerate_permutations(dataset, smiles):
    # RDKit pre-2018 used "[*]"; this changed to using a bare "*".
    if "[*]" in smiles:
        wildcard_pat = _bracket_wildcard_pat
        wildcard = "[*]"
    elif "*" in smiles:
        wildcard_pat = _organic_wildcard_pat
        wildcard = "*"
        
    n = smiles.count("*")
    if n == 1:
        yield "1", smiles.replace(wildcard, "[*:1]")
        return
            
    if n == 2:
        sub_terms = ["[*:2]", "[*:1]"]
        yield "12", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        if dataset.is_symmetric:
            return
        sub_terms = ["[*:1]", "[*:2]"]
        yield "21", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        return
    
    if n == 3:
        sub_terms = ["[*:3]", "[*:2]", "[*:1]"]
        yield "123", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        if dataset.is_symmetric:
            return
        
        sub_terms = ["[*:2]", "[*:3]", "[*:1]"]
        yield "132", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        
        sub_terms = ["[*:3]", "[*:1]", "[*:2]"]
        yield "213", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        
        sub_terms = ["[*:1]", "[*:3]", "[*:2]"]
        yield "231", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        
        sub_terms = ["[*:2]", "[*:1]", "[*:3]"]
        yield "312", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        
        sub_terms = ["[*:1]", "[*:2]", "[*:3]"]
        yield "321", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        
        return

    raise AssertionError(smiles)

# The LHS only has "*", the RHS has "*:1", "*:2", ...

_weld_cache = {}
def weld_fragments(frag1, frag2):
    key = (frag1, frag2)
    value = _weld_cache.get(key, None)
    if value is not None:
        return value

    # Also cache the lhs and rhs parts because they can be reused.
    # (It's about 4% faster overall runtime on one test.)
    
    frag1_closures = _weld_cache.get(frag1, None)
    if frag1_closures is None:
        frag1_closures = smiles_syntax.convert_wildcards_to_closures(frag1, [1, 2, 3])
        _weld_cache[frag1] = frag1_closures
        
    frag2_closures = _weld_cache.get(frag2, None)
    if frag2_closures is None:
        frag2_closures = smiles_syntax.convert_labeled_wildcards_to_closures(frag2)
        _weld_cache[frag2] = frag2_closures

    welded_mol = Chem.MolFromSmiles(frag1_closures + "." + frag2_closures)
    assert welded_mol is not None, (frag1, frag2, frag1_closures + "." + frag2_closures)
    welded_smiles = Chem.MolToSmiles(welded_mol, isomericSmiles=True)
    
    if len(_weld_cache) > 3000:
        _weld_cache.clear()
        _weld_cache[frag1] = frag1_closures
        _weld_cache[frag2] = frag2_closures
    value = (welded_smiles, welded_mol)
    _weld_cache[key] = value
    return value


def _weld_and_filter(item):
    frag_constant_smiles, frag_variable_smiles, substructure_pat, row = item
    rule_id, rule_environment_id, other_constant_smiles, is_reversed = row
    product_smiles, new_mol = weld_fragments(frag_constant_smiles, str(other_constant_smiles))
    if substructure_pat is not None:
        # The input SMARTS can contain an explict [H],
        # which in SMARTS only matches explicit hydrogens,
        # not implicit hydrogens. It's easier to make all
        # of the hydrogens explicit than it is to adjust
        # any explicit [H] terms in the query.
        test_mol = Chem.AddHs(new_mol)
        passed_substructure_test = test_mol.HasSubstructMatch(substructure_pat)
    else:
        passed_substructure_test = True
    return (frag_constant_smiles, frag_variable_smiles, row, product_smiles, passed_substructure_test)


def make_transform(
        dataset, transform_fragments, property_info_list,
        rule_selection_function,
        substructure_pat=None,
        min_radius=0, min_pairs=0, min_variable_size=0, min_constant_size=0,
        pool=None,
        cursor=None, explain=None):
    if explain is None:
        explain = command_support.no_explain
    if cursor is None:
        cursor = dataset.get_cursor()
    assert min_radius in (0, 1, 2, 3, 4, 5)

    # Map from the destination SMILES to the set of rule environments
    # The RHS set contains (rule_id, rule_environment_id, is_reversed) tuples.
    product_rule_environment_table = defaultdict(set)
    
    # Hold the welded molecules in case I need them for a substructure search
    
    # For each variable fragment (single, double, or triple cut) and
    # for each environment, extract all rules from the DB that start
    # with the given fragment and that has the same environment as the
    # query fragment (for all environment levels). Consider only
    # environments with radius >= the radius given as input argument.

    to_weld = []
    
    # This includes the possible fragments of hydrogens
    for frag in transform_fragments:
        ## Note on terminology:
        # constant = [*:1]Br.[*:2]c1ccccc1
        # variable = c1ccc(-c2sc(-c3ccc([*:1])cc3)pc2[*:2])cc1

        explain("Processing fragment %r", frag)

        # Check if the fragmentation is allowed
        if min_variable_size and frag.num_variable_heavies < min_variable_size:
            explain("  The %d heavy atoms of constant %r is below the --min-constant-size of %d. Skipping fragment.",
                    frag.num_variable_heavies, frag.variable_smiles, min_variable_size)
            continue
        if min_constant_size and frag.num_constant_heavies < min_constant_size:
            explain("  The %d heavy atoms of constant %r is below the --min-constant-size of %d. Skipping fragment.",
                    frag.num_constant_heavies, frag.constant_smiles, min_constant_size)
            continue
            
        # XXX TODO: handle 'constant_with_H_smiles'?

        # The variable SMILES contains unlabeled attachment points, while the
        # rule_smiles in the database contains labeled attachment points.
        # The fragment [*]CO[*] can potentially match [*:1]CO[*:2] or [*:2]CO[*:1],
        # so I need to enumerate all n! possibilities and find possible matches.
        
        query_possibilities = []
        for permutation, permuted_variable_smiles in enumerate_permutations(dataset, frag.variable_smiles):
            permuted_variable_smiles_id = dataset.get_rule_smiles_id(permuted_variable_smiles, cursor=cursor)
            if permuted_variable_smiles_id is not None:
                explain("  variable %r matches SMILES %r (id %d)", 
                        frag.variable_smiles, permuted_variable_smiles, permuted_variable_smiles_id)
                query_possibilities.append( (permutation, permuted_variable_smiles, permuted_variable_smiles_id) )
            else:
                explain("  variable %r not found as SMILES %r",
                        frag.variable_smiles, permuted_variable_smiles)

        if not query_possibilities:
            explain("  No matching rule SMILES found. Skipping fragment.")
            continue

        explain(" Evaluating %d possible rule SMILES: %s", 
                len(query_possibilities), sorted(x[0] for x in query_possibilities))

        # We now have a canonical variable part, and the assignment to the constant part.
        # Get its fingerprints.
        
        all_center_fps = environment.compute_constant_center_fingerprints(
            frag.constant_smiles, min_radius=min_radius)

        # For each possible way to represent the variable SMILES:
        #   Find all of the pairs which use the same SMILES id as the variable
        #   (The pairs are ordered so the matching SMILES is the 'from' side of the transform)
        #   The transformed SMILES goes from variable+constant -> dest_smiles+constant
        #   so weld the destinition SMILES (smi2) with the constant
        
        for permutation, permuted_variable_smiles, permuted_variable_smiles_id in query_possibilities:
            explain(" Evaluate constant %r with permutation %r against rules using SMILES %s (%d)",
                    frag.constant_smiles, permutation, permuted_variable_smiles, permuted_variable_smiles_id)

            possible_envs = environment.get_all_possible_fingerprints(
                all_center_fps, frag.variable_symmetry_class, permutation)
            
            rows = dataset.find_rule_environments_for_transform(
                permuted_variable_smiles_id, sorted(possible_envs), cursor=cursor)

            to_weld.extend( (frag.constant_smiles, frag.variable_smiles, substructure_pat, row) for row in rows )
        
    if pool is None:
        results = _compat.imap(_weld_and_filter, to_weld)
    else:
        # A chunk size of 20 seems to maximimize performance.
        # Too small and there's extra pickling overhead. (Larger chunks share the same SMARTS pickle.)
        # Too large and only one process might be used for all of the welding.
        results = pool.imap(_weld_and_filter, to_weld, 20)
        
    for frag_constant_smiles, frag_variable_smiles, row, product_smiles, passed_substructure_test in results:
        rule_id, rule_environment_id, other_constant_smiles, is_reversed = row
        if not passed_substructure_test:
            explain("     Skip rule %d:  %r + %r -> %r; does not contain --substructure",
                    rule_id, frag_constant_smiles, str(other_constant_smiles), product_smiles)
            continue
            
            # How to get to product_smiles from variable_smiles using rule_environment_id
        product_rule_environment_table[product_smiles].add(
            (rule_id, frag_variable_smiles, rule_environment_id, is_reversed))
        explain("     Rule %d:  %r + %r -> %r",
                rule_id, frag_constant_smiles, str(other_constant_smiles), product_smiles)
                

    explain("== Product SMILES in database: %d ==" % (len(product_rule_environment_table),))

    transform_products = list(iter_transform_products(
        dataset, product_rule_environment_table, property_info_list,
        min_pairs, rule_selection_function,
        cursor, explain))
    
    return TransformResult(property_info_list, transform_products)

class TransformProduct(object):
    def __init__(self, smiles, property_rules):
        self.smiles = smiles
        self.property_rules = property_rules

def iter_transform_products(
        dataset, product_rule_environment_table, property_info_list,
        min_pairs, rule_selection_function,
        cursor, explain):
    # For each SMILES (row of output)
    for product_smiles, rule_environment_info in sorted(product_rule_environment_table.items()):
        explain("Evaluating the %d rule environments for %r", len(rule_environment_info), product_smiles)
        product_property_rules = []

        # For each property...
        for property_name_id, property_name in property_info_list:
            explain("  Evaluating property %s (%d) for SMILES %s",
                    property_name, property_name_id, product_smiles)

            # Figure out the rule environments
            property_rules = []
            for (rule_id, variable_smiles, rule_environment_id, is_reversed) in rule_environment_info:
                property_rule = dataset.get_property_rule(
                    property_name_id, rule_environment_id, is_reversed, cursor=cursor)
                if property_rule is None:
                    explain("    No property values found for %s environment (id %d)",
                            property_name, rule_environment_id)
                    continue

                if property_rule.count < min_pairs:
                    explain("    Rule %s (#%d) for %s environment (id %d) count is too small: %d < %d",
                            property_rule.smirks, property_rule.rule_id, property_name,
                            rule_environment_id, property_rule.count, min_pairs)
                    continue
                property_rules.append(property_rule)

            if not property_rules:
                explain("    No rules to select.")
                property_rule = None
            else:
                property_rule = rule_selection_function(property_rules, explain)
                if property_rule is None:
                    explain("    No rule selected.")
                    
            product_property_rules.append(property_rule)

        yield TransformProduct(product_smiles, product_property_rules)

class TransformResult(object):
    missing_value = ""
    field_names = ("from_smiles", "to_smiles", "radius", "fingerprint", "count", "avg", "std", "kurtosis", "skewness",
                   "min", "q1", "median", "q3", "max", "paired_t", "p_value")
    
    column_formatters = _global_column_formatters.copy()
    
    def __init__(self, property_info_list, transform_products):
        self.property_info_list = property_info_list
        self.transform_products = transform_products
        
    def _write_header(self, product_file, column_names):
        # Output column names
        #   ID SMILES transformation_SMILES count_MP avg_MP ...
        output_column_names = ["ID", "SMILES"]
        for id, property_name in self.property_info_list:
            for column_name in column_names:
                output_column_names.append("%s_%s" % (property_name, column_name))
        product_file.write("\t".join(output_column_names) + "\n")
       
    def write_products(self, product_file, field_names=None,
                       column_aliases=None, include_empty=False):
        if field_names is None:
            field_names = TransformResult.field_names
        column_names = _get_column_names(field_names, column_aliases)
        
        self._write_header(product_file, column_names)

        product_id = 1
        for transform_product in self.transform_products:
            values = [None, transform_product.smiles]
            has_nonempty = False
            for property_rule in transform_product.property_rules:
                if property_rule is None:
                    rule_values = [self.missing_value]*len(field_names)
                else:
                    rule_values = _format_object_attributes(
                        property_rule, field_names, self.missing_value, self.column_formatters)
                    has_nonempty = True
                values.extend(rule_values)
            
            if has_nonempty or include_empty:
                values[0] = str(product_id)
                product_id += 1
                product_file.write("\t".join(values) + "\n")

def test_transform():
    import sys
    from . import dbutils
    db = dbutils.open_database("e.mmpdb")
    dataset = db.get_dataset()
    transform_tool = get_transform_tool(dataset)
    transform_record = transform_tool.fragment_transform_smiles("c1ccccc1C(=O)N(C)C")
    result = transform_tool.transform(transform_record.fragments, ["MP"],
                                      explain=command_support.get_explain(False))
    result.write_products(sys.stdout, include_empty=True)
    
    
if __name__ == "__main__":
    test_predict()
    test_transform()
