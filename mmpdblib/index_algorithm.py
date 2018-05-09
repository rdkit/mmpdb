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

from __future__ import absolute_import, print_function, division

from collections import defaultdict, OrderedDict
from scipy import stats
import numpy as np
import re
import sys
import os
import itertools
import json
import binascii
import operator

from . import fragment_io
from . import fileio, reporters
from . import environment
from . import fragment_algorithm
from . import index_writers
from . import _compat

###

MAX_RADIUS = 5 # maximum allowed environment radius

def _positive_float(value):
    value = float(value)
    if value <= 0.0:
        raise ValueError("must be a positive float")
    return value

def _nonnegative_float(value):
    value = float(value)
    if value < 0.0:
        raise ValueError("must be a positive float or zero")
    return value

def _nonnegative_int(value):
    value = int(value)
    if not (value >= 0):
        raise ValueError("must be a positive integer or zero")
    return value

parse_min_variable_heavies_value = _nonnegative_int
parse_max_variable_heavies_value = _nonnegative_int

parse_max_variable_ratio_value = _nonnegative_float
parse_min_variable_ratio_value = _positive_float
parse_max_heavies_transf = _nonnegative_int

class IndexOptions(object):
    def __init__(self, min_variable_heavies=None, max_variable_heavies=None,
                     min_variable_ratio=None, max_variable_ratio=None,
                     symmetric=False, max_heavies_transf=None, max_frac_trans=None):
        self.min_variable_heavies = min_variable_heavies
        self.max_variable_heavies = max_variable_heavies
        self.min_variable_ratio = min_variable_ratio
        self.max_variable_ratio = max_variable_ratio
        
        self.symmetric = symmetric
        self.max_heavies_transf = max_heavies_transf
        self.max_frac_trans = max_frac_trans
    
    def to_dict(self):
        from collections import OrderedDict
        items = [
            ("min_variable_heavies", self.min_variable_heavies),
            ("max_variable_heavies", self.max_variable_heavies),
            ("min_variable_ratio", self.min_variable_ratio),
            ("max_variable_ratio", self.max_variable_ratio),
            ("symmetric", self.symmetric),
            ("max_heavies_transf", self.max_heavies_transf),
            ("max_frac_trans", self.max_frac_trans)]
        return OrderedDict([(key, value) for (key, value) in items if value is not None])

### Filter the fragments coming in


class MinVariableRatioFilter(object):
    def __init__(self, ratio):
        self.ratio = ratio

    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        variable_ratio = num_variable_heavies / num_normalized_heavies
        return variable_ratio >= self.ratio

    def get_args(self):
        return {"--min-variable-ratio": str(self.ratio)}

    def get_options(self):
        return {"min_variable_ratio": self.ratio}

class MaxVariableRatioFilter(object):
    def __init__(self, ratio):
        self.ratio = ratio

    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        variable_ratio = num_variable_heavies / num_normalized_heavies
        return variable_ratio <= self.ratio

    def get_args(self):
        return {"--max-variable-ratio": str(self.ratio)}
        
    def get_options(self):
        return {"max_variable_ratio": self.ratio}

class MinVariableHeaviesFilter(object):
    def __init__(self, min_size):
        self.min_size = min_size
        
    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        return num_variable_heavies >= self.min_size

    def get_args(self):
        return {"--min-variable-heavies": str(self.min_size)}

    def get_options(self):
        return {"min_variable_heavies": self.min_size}

class MaxVariableHeaviesFilter(object):
    def __init__(self, max_size):
        self.max_size = max_size
        
    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        return num_variable_heavies <= self.max_size

    def get_args(self):
        return {"--max-variable-heavies": str(self.max_size)}
        
    def get_options(self):
        return {"max_variable_heavies": self.max_size}

            
class MultipleFilters(object):
    # Boolean 'and' of all of the tests
    def __init__(self, filters):
        self.filters = filters

    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        for filter in self.filters:
            if not filter.allow_fragment(num_variable_heavies, num_normalized_heavies):
                return False
        return True
    
    def get_args(self):
        d = {}
        for filter in self.filters:
            d.update(filter.get_args())
        return d

    def get_options(self):
        d = {}
        for filter in self.filters:
            d.update(filter.get_options())
        return d
    
class _AllowAllFilter(object):
    def allow_fragment(self, num_variable_heavies, num_normalized_heavies):
        return True
    

#########

# Input:
#   Two canonical SMILES with unlabeled 
# Canonical mapping from smiles1 -> smiles2
#
def sym11(p):
    return [(p[0], p[1]), (p[1], p[0])]

def sym12(p):
    return [(p[0], p[1])]

def sym111(p):
    return [(p[0], p[1], p[2]),
            (p[0], p[2], p[1]),
            (p[1], p[0], p[2]),
            (p[1], p[2], p[0]),
            (p[2], p[1], p[0]),
            (p[2], p[0], p[1])]
def sym112(p):
    return [(p[0], p[1], p[2]),
            (p[1], p[0], p[2])]
def sym121(p):
    return [(p[0], p[1], p[2]),
            (p[2], p[1], p[0])]
def sym122(p):
    return [(p[0], p[1], p[2]),
            (p[0], p[2], p[1])]
def sym123(p):
    return [(p[0], p[1], p[2])]

_symm_funcs = {
    "11": sym11,
    "12": sym12,
    "111": sym111,
    "112": sym112,
    "121": sym121,
    "122": sym122,
    "123": sym123,
    }
def _check_sym_funcs():
    seen_results = set()
    

    for name, f in _symm_funcs.items():
        n = len(name)
        query = list(range(4, 4+n))
        result = f(query)
        
        x = tuple(sorted(result))
        assert x not in seen_results, x
        seen_results.add(x)
        
        seen_terms = set()
        for term in result:
           assert len(term) == n, (query, term)
           assert set(term) == set(query), (query, term)
           assert term not in seen_terms
           seen_terms.add(term)
_check_sym_funcs()        

def enumerate_symmetry(possibilities, symmetry_class):
    func = _symm_funcs[symmetry_class]
    result = []
    for p in possibilities:
        for q in func(p[-1]):
            new_p = p + (q,)
            result.append(new_p)
    return result

def reorder(possibilities, order):
    result = []
    for p in possibilities:
        q = tuple(p[-1][offset] for offset in map(int, order))
        new_p = p + (q,)
        result.append(new_p)
    return result

def _invert_order(order):
    values = list(map(int, order))
    indices = list(range(len(values)))
    indices.sort(key = lambda i: values[i])
    return "".join(str(i) for i in indices)


_order_table = dict((s, _invert_order(s))
                       for s in ("01", "10", "012", "021", "102", "120", "201", "210"))

def invert_order(order):
    return _order_table[order]    

def _get_smirks_order(
        symmetry_class1,
        attachment_order1,
        constant_symmetry_class,
        symmetry_class2,
        attachment_order2):
    # the i-th variable * goes to the attachment_order[i]-th constant *.
    # "120" will generate ...[*:2]...[*:3]...[*:1]... for the variable part.
    # The constant part will always be :1, :2, :3

    # That attachment order is canonical with respect to the constant part.
    # However, different constants involving the same variable (smirks)
    # part may generate different attachment orders. Consider:
    #  input   variable  constant  attachment order
    #   CXYB     *XY*      *B.*C     10
    #   BXXC     *XX*      *B.*C     01
    #   CXYD     *XY*      *C.*D     01
    #   CXXD     *XX*      *C.*D     01
    #
    # A simple canonicalization of the first pair maps
    #    *XY* to [*:2]XY[*:1] and *XX* to [*:1]XX[*:2]
    #    giving [*:2]XY[*:1]>>[*:1]XX[*:2]
    # while a canonicalization of the second pair maps
    #    *XY* to [*:1]XY[*:2] and *XX* to [*:1]XX[*:2]
    #    giving [*:1]XY[*:2]>>[*:1]XX[*:2]
    #
    # That's why the simple mapping doesn't work.
    # Instead, always use the order [*:1], [*:2], [*:3]
    # for the LHS and map them to the lowest-ordered
    # equivalent wildcard on the RHS.

    indices = tuple(range(len(symmetry_class1)))
    possibilities = [(indices,)]
    #print("possibilities0", possibilities)
    #print("symmetry_class1:", symmetry_class1, "symmetry_class2:", symmetry_class2)
    possibilities = enumerate_symmetry(possibilities, symmetry_class1) # symmetry in smiles1
    #print("possibilities1", possibilities)
    possibilities = reorder(possibilities, invert_order(attachment_order1))  # map to the constant
    #print("possibilities2", possibilities)
    possibilities = enumerate_symmetry(possibilities, constant_symmetry_class)  # symmetry in the constant
    #print("possibilities3", possibilities)
    possibilities = reorder(possibilities, attachment_order2)              # map to smiles2
    #print("possibilities4", possibilities)
    possibilities = enumerate_symmetry(possibilities, symmetry_class2)     # symmetry in smiles2
    #print("possibilities5", possibilities)
    
    #print(sorted(possibilities))
    possibilities = [(p[5], p[3]) for p in possibilities]
    
    best_possibility = min(possibilities)
    
    s = "".join(str(i) for i in best_possibility[0])
    t = "".join(str(i) for i in best_possibility[1])
    return s, t

def _init_cansmirks_table():
    # Enumerate all possibilities. Use brute-force to minimize errors.

    # There will be 4532 elements.
    table = {}
    for symmetry_classes, permutations in (
            (("11", "12"), ("01", "10")),
            (("111", "112", "121", "122", "123"),
             ("012", "021", "102", "120", "201", "210")),
             ):

        for symmetry_class1 in symmetry_classes:
            for attachment_order1 in permutations:
                for constant_symmetry_class in symmetry_classes:
                    for attachment_order2 in permutations:
                        for symmetry_class2 in symmetry_classes:
                            smirks_order = _get_smirks_order(
                                symmetry_class1,
                                attachment_order1,
                                constant_symmetry_class,
                                symmetry_class2,
                                attachment_order2)
                            
                            table[symmetry_class1 +
                                  attachment_order1 +
                                  constant_symmetry_class +
                                  symmetry_class2 +
                                  attachment_order2] = smirks_order
    #print(len(table), "elements")
    return table

# Using the pre-computed table is almost 3x faster, but it's
# harder to debug failures in smirks canonicalization.
USE_SMIRKS_TABLE = True
USE_PRECOMPUTED_TABLE = True

if USE_SMIRKS_TABLE:
    if USE_PRECOMPUTED_TABLE:
        # Saves 1.5 seconds on my computer.
        from .cansmirks_table import cansmirks_table as _smirks_table
    else:
        _smirks_table = _init_cansmirks_table()
        if 0:
            # Code used to generated the pre-computed table.
            import pprint
            with open("cansmirks_table.py", "w") as outfile:
                outfile.write("# The contents of this file were generated from index_algorithm.py. Do not modify.\n")
                outfile.write("#pylint: disable=bad-continuation\n")
                outfile.write("cansmirks_table = ")
                pprint.pprint(_smirks_table, stream=outfile)
            raise SystemExit("cansmirks_table.py generated")
    

class RelabelCache(dict):
    def __missing__(self, key):
        if isinstance(key, _compat.basestring):
            result = fragment_io.relabel(key)
        else:
            smiles, order = key
            result = fragment_io.relabel(smiles, order)
        self[key] = result
        return result

_wildcard_regex = re.compile(re.escape("[*]") + "|" + re.escape("*"))

def cansmirks(num_cuts,
              smiles1, symmetry_class1, attachment_order1,
              constant_smiles, constant_symmetry_class,
              smiles2, symmetry_class2, attachment_order2,
              relabel_cache
              ):
    if num_cuts == 1:
        # This is easy enough that I'll relabel them directly
        smirks = smiles1 + ">>" + smiles2
        if "*" in smirks:
            smirks = _wildcard_regex.sub(lambda x: "[*:1]", smirks)
        else:
            raise AssertionError(smirks)
        if "*" in constant_smiles:
            constant_smiles = _wildcard_regex.sub(lambda x: "[*:1]", constant_smiles)
        else:
            raise AssertionError(constant_smiles)
        return smirks, constant_smiles

    if USE_SMIRKS_TABLE:
        new_order, constant_order = _smirks_table[
            symmetry_class1 + attachment_order1 +
            constant_symmetry_class +
            symmetry_class2 + attachment_order2]
    else:
        ## print("smiles1:", smiles1, "smiles2:", smiles2)
        new_order, constant_order = _get_smirks_order(
                                    symmetry_class1,
                                    attachment_order1,
                                    constant_symmetry_class,
                                    symmetry_class2,
                                    attachment_order2)
        ## print("new_order:", new_order, "constant_order:", constant_order)

    ## print("Relabel", constant_smiles, "with", constant_order)
    ## print("smiles1", smiles1, " -> ", fragment_io.relabel(smiles1))
    ## print("smiles2", smiles2, " -> ", fragment_io.relabel(smiles2, new_order), "new_order:", new_order)
    ## print("constant", fragment_io.relabel(constant_smiles, constant_order))
    
    return (relabel_cache[smiles1] + ">>" + relabel_cache[smiles2, new_order],
                relabel_cache[constant_smiles, constant_order])
    

class FragmentIndex(object):
    def __init__(self, index, id_to_record):
        self._index = index
        self._id_to_record = id_to_record

    def __len__(self):
        n = sum(1 for matches in self._index.values() if len(matches) > 1)
        return n
        
    def iter_constant_matches(self):
        for (constant_smiles, constant_symmetry_class, num_cuts), matches in sorted(self._index.items()):
            if len(matches) <= 1:
                continue
            yield num_cuts, constant_smiles, constant_symmetry_class, matches

    def get_input_record(self, id):
        return self._id_to_record[id]

        

class InputRecord(object):
    __slots__ = ("id", "input_smiles", "num_normalized_heavies", "normalized_smiles")
    def __init__(self, id, input_smiles, num_normalized_heavies, normalized_smiles):
        self.id = id
        self.input_smiles = input_smiles
        self.num_normalized_heavies = num_normalized_heavies
        self.normalized_smiles = normalized_smiles
    

def load_fragment_index(fragment_reader, fragment_filter=None, selected_ids=None):
    if fragment_filter is None:
        fragment_filter = _AllowAllFilter()

    # constant -> list of records containing that constant
    index = defaultdict(list)

    # SMILES for the constant (the constant fragment R-groups) -> all corresponding identifiers
    # (There may be multiple identifiers due to de-salting.)
    normalized_smiles_to_ids = defaultdict(list)
    id_to_record = {}

    constant_smiles_to_hydrogen_constant_smiles = {}
    
    for recno, record in enumerate(fragment_reader, 1):
        if record.errmsg:
            continue
        if record.id in id_to_record:
            raise ValueError("Duplicate identifier %r at %s" % (record.id, fragment_reader.location.where()))
        id_to_record[record.id] = InputRecord(
            record.id, record.input_smiles, record.num_normalized_heavies, record.normalized_smiles)
        if selected_ids is not None and record.id not in selected_ids:
            continue
        
        normalized_smiles_to_ids[record.normalized_smiles].append(record.id)
        
        for fragmentation in record.fragments:
            if not fragment_filter.allow_fragment(
                    fragmentation.variable_num_heavies, record.num_normalized_heavies):
                continue

            index[(fragmentation.constant_smiles, fragmentation.constant_symmetry_class, fragmentation.num_cuts)].append(
                (record.id, fragmentation.variable_symmetry_class,
                 fragmentation.variable_smiles, fragmentation.attachment_order,
                 fragmentation.enumeration_label))
            
            if fragmentation.num_cuts == 1:
                constant_smiles_to_hydrogen_constant_smiles[
                    fragmentation.constant_smiles] = fragmentation.constant_with_H_smiles

    ## Add the single cut hydrogen transformations

    # The algorithm is:
    #   - for each single cut constant, get its with-hydrogen version
    #   - if the with-hydrogen version matches an actual record
    #   - add the records using the [*:1][H] variable fragment
    
    for (constant_smiles, constant_symmetry_class, num_cuts), matches in index.items():
        if num_cuts != 1:
            continue
        constant_with_H_smiles = constant_smiles_to_hydrogen_constant_smiles.get(constant_smiles, None)
        if constant_with_H_smiles is None:
            continue
        other_ids = normalized_smiles_to_ids.get(constant_with_H_smiles, [])
        for other_id in other_ids:
            # NOTE: this is hard-coded to "[*][H]", and must match the
            # same string used in fragment_algorithm.py's _hydrogen_cut_smiles
            matches.append( (other_id, "1", "[*][H]", "0", "N") )


    return FragmentIndex(dict(index), id_to_record)


class MatchedMolecularPair(object):
    __slots__ = ("id1", "id2", "smirks", "constant_smiles", "max_constant_radius")
    def __init__(self, id1, id2, smirks, constant_smiles, max_constant_radius):
        self.id1 = id1
        self.id2 = id2
        self.smirks = smirks
        self.constant_smiles = constant_smiles
        self.max_constant_radius = max_constant_radius
    

## I needed a way to get the number of heavy atoms in the variable
## fragment so I could implement the --max-heavies-transf option. This
## was done very late in the process. The easiest solution is to use a
## regex to count the number of heavy atoms in the SMILES.  The
## cleaner solution would be to pass that information all the way
## through the system.
# Extract just the atom terms (no closures)
_atom_pat = re.compile(r"""
(
 Cl? |
 Br? |
 [NOSPFIbcnosp] |
 \[[^]]*\]
)
""", re.X)

def get_num_heavies(smiles):
    num_atoms = 0
    for m in _atom_pat.finditer(smiles):
        # I think a cleverer pattern could exclude these matches.
        # \[  (\d+[A-Z][^]]*\]) | ([AB..GI...Y][^]]*\]) | (H[^]]+\])
        #     ^^^^^ has isotope for cases like [2H], and not a *
        #                     ^^^^^ element which isn't H or *
        #                         odd cases like [H+]      ^^^^^^
        # TODO: evaluate the performance
        text = m.group()
        if text == "[H]" or "*" in text:
            continue
        num_atoms += 1
    return num_atoms

assert 1/2 != 0, "why did you disable future division?"

# A different approach to caching
class _NumHeaviesCache(dict):
    def __missing__(self, smiles):
        n = get_num_heavies(smiles)
        self[smiles] = n
        return n
_num_heavies_cache = _NumHeaviesCache()

def get_max_radius_for_fraction_transfer(
        max_frac_trans, smirks, constant_smiles):
    # Need to figure out the maximum radius to use

    # A pair has compounds C1 and C2. There are N1 and N2 atoms in each, respectively.
    # The constant part has 'n' atoms.
    #     => The variable part of C1 has V1=N1-n atoms, and the variable part of C2 has V2=N2-n atoms
    #     => V1, V2 are the heavies on the LHS and RHS of the smirks

    # The constant part has circular fingerprints for R0, R1, R3, R3 ...,
    # with M0, M1, M2, M3, ... heavy atoms. (By construction, M0 will always be 0 because
    # it's the circular fingerprints of 1, 2, or 3 "*" atoms.)

    # Need to find the largest r such that
    #   max((V1+Mr)/N1, (V2+Mr)/N2) <= max_frac_trans
    #  -or-
    #   max((N1-n+Mr)/N1, (N2-n+Mr)/N2) <= max_frac_trans

    # Trivial case
    if max_frac_trans >= 1.0:
        return 0
    variable_smiles1, _, variable_smiles2 = smirks.partition(">>")
    assert _ == ">>", smirks
    
    V1 = _num_heavies_cache[variable_smiles1]
    V2 = _num_heavies_cache[variable_smiles2]
    n = _num_heavies_cache[constant_smiles]
    N1 = V1 + n
    N2 = V2 + n

    if N1 == 0 or N2 == 0:
        # Might happen if the variable is a hydrogen?
        # Don't want an divide-by-zero error if something odd happens
        return 0

    # Is r=0 allowed?
    if max(V1/N1, V2/N2) > max_frac_trans:
        return None
    
    # Quick reject test for r=1 is even allowed. There must be at least 1 heavy.
    # (It's hard to be more clever. The constant SMILES might be [*]C([*])CCC,
    # where r=1 has only one heavy atom.)
    
    if max((V1-1)/N1, (V2-1)/N2) > max_frac_trans:
        return 0

    # Otherwise we need to compute the number of atoms in the circular regions
    # XXX This is an unlimited cache, which shares global state!
    # With max_frac_trans = 1.0000, t = 0:51 
    # With max_frac_trans = 0.9999, t = 1:58 using get_cached_centers()
    # With max_frac_trans = 0.9999, t = 1:07 using get_cached_center_radii()
    #    Nearly all of that extra time is in get_num_heavies()! (Measured it at 16 seconds!)
    # With max_frac_trans = 0.9999, t = 0:53 using get_cached_center_radii() and _NumHeaviesCache()
    #    Only a few seconds longer. I think that's good enough.
    
    # If I just cache the center then my test case took 
    radii = environment_cache.get_or_compute_center_radii(constant_smiles, MAX_RADIUS)

    ## print("get_max_radius_for_fraction_transfer()")
    ## print("n", n, "N1", N1, "V1", V1, "N2", N2, "V2", V2)
    best_radius = None
    for radius, num_atom_in_radius in radii:
        frac_trans1 = (V1+num_atom_in_radius)/N1
        frac_trans2 = (V2+num_atom_in_radius)/N2
        ## print("radius", radius, "num_atom_in_radius", num_atom_in_radius,
        ##           "frac_trans1", frac_trans1, "frac_trans2", frac_trans2)
        frac_trans = min(frac_trans1, frac_trans2)
        if frac_trans <= max_frac_trans:
            best_radius = radius
        else:
            break

    ## print("best", smirks, constant_smiles, best_radius)
    return best_radius



class EnvironmentCache(object):
    def __init_(self, index_cache):
        self.index_cache = index_cache
        
    def __init__(self):
        self._centers_cache = {}
        self._radii_cache = {}
        self._environment_cache = {}
        self._interned_fingerprints = {}

    def get_or_compute_centers(self, constant_smiles):
        centers = self._centers_cache.get(constant_smiles, None)
        if centers is None:
            centers = environment.find_centers(constant_smiles.encode("ascii"))
            self._centers_cache[constant_smiles] = centers
        return centers

    def get_or_compute_center_radii(self, constant_smiles, max_radius):
        key = (constant_smiles, max_radius)
        radii = self._radii_cache.get(key, None)
        if radii is None:
            centers = self.get_or_compute_centers(constant_smiles)
            radii = list(enumerate(environment.iter_num_atoms_for_radii(centers, max_radius)))
            self._radii_cache[key] = radii
        return radii

    def get_or_compute_constant_environment(self, constant_smiles, max_radius):
        key = (constant_smiles, max_radius)
        env_fps = self._environment_cache.get(key, None)
        if env_fps is None:
            centers = self.get_or_compute_centers(constant_smiles)
            env_fps = environment.compute_constant_environment_from_centers(centers, max_radius)
            assert len(env_fps) == (max_radius+1), (len(env_fps), max_radius)
            # Many fingerprints are duplicates. Use an intern dictionary to
            # reduce the memory used. (Is this really needed/useful?)
            for env_fp in env_fps:
                fp = env_fp.fingerprint
                env_fp.fingerprint = self._interned_fingerprints.setdefault(fp, fp)
                
            self._environment_cache[key] = env_fps
        return env_fps

    

def find_matched_molecular_pairs(
        index, index_options=IndexOptions(), 
        reporter=None):

    symmetric = index_options.symmetric
    max_heavies_transf = index_options.max_heavies_transf
    max_frac_trans = index_options.max_frac_trans
    
    counter = itertools.count(0)
    reporter = reporters.get_reporter(reporter)

    relabel_cache = RelabelCache()
    NO_ENUMERATION = fragment_algorithm.EnumerationLabel.NO_ENUMERATION

    with reporter.progress(
            index.iter_constant_matches(), "Constant fragment matches", len(index)) as it:
        # Go through the upper-diagonal matrix of the NxN matches
        for num_cuts, constant_smiles, constant_symmetry_class, matches in it:
            for offset, (id1, symmetry_class1, smiles1, attachment_order1, enumeration_label1) in enumerate(matches):
                for (id2, symmetry_class2, smiles2, attachment_order2, enumeration_label2) in matches[offset+1:]:
                    ## print("QQQ")
                    ## print("1:", id1, symmetry_class1, smiles1, attachment_order1, enumeration_label1)
                    ## print("2:", id2, symmetry_class2, smiles2, attachment_order2, enumeration_label2)
                    # Eliminate duplicates based on id or matching variable fragment
                    if id1 == id2:
                        continue

                    if max_heavies_transf is not None:
                        num_heavies1 = get_num_heavies(smiles1)
                        num_heavies2 = get_num_heavies(smiles2)
                        if abs(num_heavies2-num_heavies1) > max_heavies_transf:
                            continue

                    # Simple rejection
                    if smiles1 == smiles2 and attachment_order1 == attachment_order2:
                        continue

                    # "Two constant parts may not be matched if both of them have the CHI_UP tag"
                    if (enumeration_label1 != NO_ENUMERATION and
                        enumeration_label2 != NO_ENUMERATION):
                        continue
                    
                    parameters = [
                        (id1, smiles1, symmetry_class1, attachment_order1,
                         id2, smiles2, symmetry_class2, attachment_order2),
                        (id2, smiles2, symmetry_class2, attachment_order2,
                         id1, smiles1, symmetry_class1, attachment_order1),
                         ]

                    # Put them in canonical order.
                    if (smiles1, attachment_order1) > (smiles2, attachment_order2):
                        parameters.reverse()

                    if not symmetric:
                        del parameters[1]

                    for (tmp_id1, tmp_smiles1, tmp_symmetry_class1, tmp_attachment_order1,
                         tmp_id2, tmp_smiles2, tmp_symmetry_class2, tmp_attachment_order2) in parameters:
                        smirks, tmp_constant_smiles = cansmirks(
                            num_cuts,
                            tmp_smiles1, tmp_symmetry_class1, tmp_attachment_order1,
                            constant_smiles, constant_symmetry_class,
                            tmp_smiles2, tmp_symmetry_class2, tmp_attachment_order2,
                            relabel_cache)

                        if 0:
                            # Double-check that the new assignments are valid
                            # This has a large overhead. For example, what took
                            # 29 seconds takes 143 seconds when this is enabled.
                            from . import smiles_syntax
                            from rdkit import Chem
                            
                            print("= Check", next(counter), "=")
                            print(" mol1:", tmp_id1, tmp_smiles1, tmp_symmetry_class1, tmp_attachment_order1)
                            print(" constant:", constant_smiles, constant_symmetry_class)
                            print(" mol2:", tmp_id2, tmp_smiles2, tmp_symmetry_class2, tmp_attachment_order2)
                            print(" smirks:", smirks)
                            print(" new constant:", tmp_constant_smiles)
                            const_smi = smiles_syntax.convert_labeled_wildcards_to_closures(tmp_constant_smiles)
                            lhs, _, rhs = smirks.partition(">>")
                            var_smi = smiles_syntax.convert_labeled_wildcards_to_closures(lhs)
                            smi = Chem.CanonSmiles(var_smi + "." + const_smi)
                            expected_smi = index.get_input_record(tmp_id1).normalized_smiles
                            print(" LHS  got:", smi, "from", var_smi+"."+const_smi)
                            print(" expected:", expected_smi)
                            if smi != expected_smi:
                                if Chem.CanonSmiles(smi, 0) != Chem.CanonSmiles(expected_smi, 0):
                                    print(" !!! FAILED BAD1 !!!", smi, expected_smi)
                                else:
                                    print(" !!! FAILED CHI1 !!!", smi, expected_smi)
                                #raise AssertionError(smi, expected_smi)
                            
                            var_smi = smiles_syntax.convert_labeled_wildcards_to_closures(rhs)
                            smi = Chem.CanonSmiles(var_smi + "." + const_smi)
                            expected_smi = index.get_input_record(tmp_id2).normalized_smiles
                            print(" RHS  got:", smi, "from", var_smi+"."+const_smi)
                            print(" expected:", expected_smi)
                            if smi != expected_smi:
                                if Chem.CanonSmiles(smi, 0) != Chem.CanonSmiles(expected_smi, 0):
                                    print(" !!! FAILED BAD2 !!!", smi, expected_smi)
                                else:
                                    print(" !!! FAILED CHI2 !!!", smi, expected_smi)
                                #raise AssertionError(smi, expected_smi)
                            print("   Good!")
                            

                        # Figure out the max radius allowed for environment fingerprints
                        if max_frac_trans is None or max_frac_trans >= 1.0:
                            max_radius = MAX_RADIUS
                        else:
                            max_radius = get_max_radius_for_fraction_transfer(
                                max_frac_trans, smirks, tmp_constant_smiles)
                            if max_radius is None:
                                # skip this pair
                                continue
                        
                        yield MatchedMolecularPair(tmp_id1, tmp_id2, smirks, tmp_constant_smiles, max_radius)

###
class BaseWriter(object):
    def __init__(self, backend, fragment_options, fragment_index, index_options, properties):
        self.backend = backend
        self.fragment_options = fragment_options
        self.fragment_index = fragment_index
        self.index_options = index_options
        self.properties = properties  # properties_io.Properties instance.
        if properties is None:
            self.num_properties = 0
        else:
            self.num_properties = len(properties.property_columns)

    def start(self):
        pass

    def end(self, reporter=None):
        pass
    
    def close(self):
        self.backend.close()
        
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        if type is not None:
            self.backend.rollback()
        else:
            self.backend.commit()

    def write_matched_molecule_pairs(self, pairs):
        raise NotImplementedError
    
    
class CSVPairWriter(BaseWriter):
    def start(self):
        self.num_pairs = 0
    
    def write_matched_molecule_pairs(self, pairs):
        backend = self.backend
        fragment_index = self.fragment_index

        n = 0
        for pair in pairs:
            rec1 = fragment_index.get_input_record(pair.id1)
            rec2 = fragment_index.get_input_record(pair.id2)
            backend.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                rec1.input_smiles, rec2.input_smiles, pair.id1, pair.id2, pair.smirks, pair.constant_smiles))
            n += 1
        self.num_pairs += n

    def __exit__(self, type, value, traceback):
        self.backend.close()

class RuleSmilesTable(dict):
    def __init__(self, backend):
        self.backend = backend

    def __missing__(self, smiles):
        idx = len(self)
        self.backend.add_rule_smiles(idx, smiles)
        self[smiles] = idx
        return idx
    
class ConstantSmilesTable(dict): # XXX Merge with SmilesTable?
    def __init__(self, backend):
        self.backend = backend

    def __missing__(self, smiles):
        idx = len(self)
        self.backend.add_constant_smiles(idx, smiles)
        self[smiles] = idx
        return idx
        
class RuleTable(dict):
    def __init__(self, rule_smiles_table, backend):
        self.rule_smiles_table = rule_smiles_table
        self.backend = backend

    def __missing__(self, smirks):
        from_smiles, gtgt, to_smiles = smirks.partition(">>")
        assert gtgt == ">>", smirks
        
        rule_idx = len(self)
        from_smiles_idx = self.rule_smiles_table[from_smiles]
        to_smiles_idx = self.rule_smiles_table[to_smiles]

        self.backend.add_rule(rule_idx, from_smiles_idx, to_smiles_idx)
        self[smirks] = rule_idx
        return rule_idx
        
        
        
class RuleEnvironmentTable(dict):
    def __init__(self, num_properties, environment_cache, backend):
        self.num_properties = num_properties
        self.environment_cache = environment_cache
        self.backend = backend

    def __missing__(self, key):
        # Note: "string-encoded-key" (see corresponding note, below)
        # Decode the key to get the actualy key data.
        rule_idx, env_fp_idx, radius = map(int, key.split(","))

        # To save space, if there are no properties then only store the index,
        # otherwise store the index and lists of property values.
        rule_env_idx = len(self)
        if self.num_properties == 0:
            rule_env = rule_env_idx
        else:
            property_value_lists = [[] for _ in range(self.num_properties)]
            rule_env = RuleEnvironment(rule_env_idx, property_value_lists)
        
        self.backend.add_rule_environment(rule_env_idx, rule_idx, env_fp_idx, radius)
        self[key] = rule_env
        return rule_env

    def iter_sorted_rule_environments(self):
        if self.num_properties == 0:
            # These are raw integers. Wrap them in an empty RuleEnvironment.
            # Note: I don't actualy use this branch in the code
            raise NotImplementedError("I haven't tested this branch because it isn't needed.")
#            indices = sorted(self.values())
#            empty_list = []
#            for rule_env_idx in indices:
#                yield RuleEnvironment(rule_env_idx, empty_list)
        else:
            rule_envs = list(self.values())
            # The operator replaces "lambda rule_env: rule_env.idx"
            rule_envs.sort(key=operator.attrgetter("idx"))
            for rule_env in rule_envs:
                yield rule_env

        
class EnvironmentFingerprintTable(dict):
    def __init__(self, backend):
        self.backend = backend

    def __missing__(self, env_fp):
        idx = len(self)
        self.backend.add_environment_fingerprint(idx, env_fp)
        self[env_fp] = idx
        return idx
    
class CompoundTable(dict):
    def __init__(self, fragment_index, property_name_idxs, properties, backend):
        self.fragment_index = fragment_index
        self.property_name_idxs = property_name_idxs
        self.properties = properties
        self.backend = backend

    def __missing__(self, compound_id):
        compound_idx = len(self)
        
        record = self.fragment_index.get_input_record(compound_id)

        # "id", "public_id", "input_smiles", "clean_smiles", "clean_num_heavies"]
        self.backend.add_compound(compound_idx, compound_id, record.input_smiles,
                                  record.normalized_smiles, record.num_normalized_heavies)
        self[compound_id] = compound_idx

        if self.properties:
            property_values = self.properties.get_property_values(compound_id)
            for property_idx, value in zip(self.property_name_idxs, property_values):
                if value is not None:
                    self.backend.add_compound_property(compound_idx, property_idx, value)
            
        return compound_idx

## _heap = None    
class MMPWriter(BaseWriter):
    def start(self):
        self._environment_pair_id_counter = itertools.count(0)
        self._environment_cache = EnvironmentCache()

        self.property_name_idxs = property_name_idxs = []
        if self.properties is not None:
            for property_name_idx, property_name in enumerate(self.properties.property_names):
                self.backend.add_property_name(property_name_idx, property_name)
                property_name_idxs.append(property_name_idx)
        
        self._rule_smiles_table = RuleSmilesTable(self.backend)
        self._constant_smiles_table = ConstantSmilesTable(self.backend)
        self._rule_table = RuleTable(self._rule_smiles_table, self.backend)
        self._fingerprint_table = EnvironmentFingerprintTable(self.backend)
        self._compound_table = CompoundTable(
            self.fragment_index, property_name_idxs, self.properties, self.backend)
        self._rule_environment_table = RuleEnvironmentTable(
            self.num_properties, self._environment_cache, self.backend)

        self.backend.start(self.fragment_options, self.index_options)
        self.num_pairs = 0
            
                
    def end(self, reporter=None):
        reporter = reporters.get_reporter(reporter)
        if self.properties is not None:
            add_rule_environment_statistics = self.backend.add_rule_environment_statistics

            for property_i, property_name_idx in enumerate(self.property_name_idxs):
                property_name = self.properties.property_names[property_i]
                with reporter.progress(self._rule_environment_table.iter_sorted_rule_environments(),
                                       "Writing rule statistics for property %s:"  % (property_name,),
                                       len(self._rule_environment_table)) as rule_env_iter:
                    for rule_env in rule_env_iter:
                        value_list = rule_env.property_value_lists[property_i]
                        if value_list:
                            add_rule_environment_statistics(
                                rule_env.idx, property_name_idx,
                                compute_aggregate_values(value_list))


        self.backend.end(reporter)

    def write_matched_molecule_pairs(self, pairs):
        # The main entry point for writing results to a file.
        
        has_properties = (self.properties is not None)

        pair_i = -1
        for pair_i, pair in enumerate(pairs):
            # Figure out which rule it goes into.
            rule_idx = self._rule_table[pair.smirks]

            # Get the environment for the constant part, at different radii.
            rule_envs = self._get_rule_environments(rule_idx, pair.constant_smiles, pair.max_constant_radius)
            if rule_envs:
                compound1_idx = self._compound_table[pair.id1]
                compound2_idx = self._compound_table[pair.id2]
                constant_idx = self._constant_smiles_table[pair.constant_smiles]

                for rule_env in rule_envs:
                    pair_idx = next(self._environment_pair_id_counter)

                    if has_properties:
                        # then the rule_env is a RuleEnvironment instance
                        self.backend.add_rule_environment_pair(
                            pair_idx, rule_env.idx, compound1_idx, compound2_idx, constant_idx)
                        rule_env.append_pair_properties(
                            self.properties.get_property_values(pair.id1),
                            self.properties.get_property_values(pair.id2))
                    else:
                        # then the rule_env is an integer
                        self.backend.add_rule_environment_pair(
                            pair_idx, rule_env, compound1_idx, compound2_idx, constant_idx)
                        
        self.num_pairs += (pair_i + 1)
        
    def _get_rule_environments(self, rule_idx, constant_smiles, max_radius):
        # XXX Add another layer of cache? I don't think it makes much sense.
        env_fps = self._environment_cache.get_or_compute_constant_environment(constant_smiles, max_radius)
        rule_envs = []
        for env_fp in env_fps:
            env_fp_idx = self._fingerprint_table[env_fp.fingerprint]
            # NOTE: "string-encoded-key"
            # Originally I stored the tuple directly. This ends up using a lot of space
            # because there can be tens of millions of rule environments. A 3-element
            # tuple needs 80 bytes on CPython 2.7 and 72 on CPython 3.6, plus the space
            # for each of the integers. The radius is small, so CPython uses a cached
            # value. Still, that's 28 or 2*28 bytes extra.
            # On the other hand, an encoded-string takes only about 70 bytes total.
            # Across 10M objects, this saves about 10E6*(110-70)/1024/1024 = 380 GB.
            key = "%d,%d,%d" % (rule_idx, env_fp_idx, env_fp.radius)
            rule_env = self._rule_environment_table[key]
            rule_envs.append(rule_env)

        return rule_envs

            
def open_mmpa_writer(destination, format, title, fragment_options,
                     fragment_index, index_options, properties,
                     environment_cache):
    if format is None:
        if destination is None:
            format = "mmpa"
        else:
            s = destination.lower()
            for suffix, format_name in (
                    (".csv", "csv"),
                    (".csv.gz", "csv.gz"),
                    (".mmpa", "mmpa"),
                    (".mmpa.gz", "mmpa.gz"),
                    (".mmpdb", "mmpdb"),
                    (".mmpz", "mmpz"), # XXX REMOVE
                    ):
                if s.endswith(suffix):
                    format = format_name
                    break
            else:
                format = "mmpdb"

    if format in ("csv", "csv.gz"):
        # Can be a helpful summary
        outfile = fileio.open_output(destination, format)
        return CSVPairWriter(outfile, fragment_options, fragment_index,
                             index_options, properties)
    elif format in ("mmpa", "mmpa.gz"):
        # Text-based version somewhat useful for debugging
        outfile = fileio.open_output(destination, format)
        index_writer = index_writers.open_table_index_writer(outfile)
    elif format == "mmpdb":
        # The preferred output format.
        index_writer = index_writers.open_sqlite_index_writer(destination, title)
        
    else:
        raise AssertionError(format)

    return MMPWriter(index_writer, fragment_options, fragment_index,
                     index_options, properties)
    
# Parameters

#OPTION key value
#COMPOUND id public_id input_smiles clean_smiles clean_num_heavies
#RULE id smirks
#PAIR id rule_id compound_id1 compound_id2

class Compound(object):
    def __init__(self, idx, compound_id, record, property_values):
        self.idx = idx
        self.compound_id = compound_id
        self.record = record
        self.property_values = property_values
        


# From https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
# with changes to return None if there isn't enough data
def online_variance(values): 
    n = 0
    mean = 0.0
    M2 = 0.0
     
    for x in values:
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    if n < 2:
        return None
    else:
        return M2 / (n - 1)

def online_kurtosis(data):
    n = 0
    mean = 0
    M2 = 0
    M3 = 0
    M4 = 0

    for x in data:
        n1 = n
        n = n + 1
        delta = x - mean
        delta_n = delta / n
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        mean = mean + delta_n
        M4 = M4 + term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3
        M3 = M3 + term1 * delta_n * (n - 2) - 3 * delta_n * M2
        M2 = M2 + term1

    if M2 == 0:
        return None
    kurtosis = (n*M4) / (M2*M2) - 3
    return kurtosis    

def get_median(values):
    n = len(values)
    if n == 0:
        return None
    
    half = n//2
    
    if n % 2 == 1:
        # Odd number, like: 5, 7, 9
        # 3//2 = 1, which is the median
        median = values[half]
    else:
        # Even number, like: 5, 7, 9, 11
        # 4//2 = 2, which is the index just above the middle
        median = (values[half-1] + values[half]) / 2
    return median

# Compute interpolated quartiles according to "Method 3" of
# https://en.wikipedia.org/wiki/Quartile#Method_3
def compute_quartiles(values):
    n = len(values)
    assert n > 0
    if n == 1:
        return (values[0], values[0], values[0])

    median = get_median(values)
    
    half = n//2
    if n % 2 == 0:
        q1 = get_median(values[:half])
        q3 = get_median(values[half:])
    elif n % 4 == 1:
        m = (n-1) // 4
        q1 = 0.25*values[m-1] + 0.75 * values[m]
        q3 = 0.75*values[3*m] + 0.25 * values[3*m+1]
    else:
        assert n % 4 == 3
        m = (n-3) // 4
        q1 = 0.75*values[m] + 0.25 * values[m+1]
        q3 = 0.25*values[3*m+1] + 0.75 * values[3*m+2]
    return q1, median, q3

if __debug__:
    for test_data, expected in (
            ([3], (3, 3, 3)),
            ([3, 4], (3, 3.5, 4)),
            ([3, 4, 5], (3.25, 4, 4.75)),
            # Test cases from Wikipedia
            ([6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49], (20.25, 40, 42.75)), # 11 = 4*2+3
            ([7, 15, 36, 39, 40, 41], (15.0, 37.5, 40.0)),  # 6 = 4*1+2
            # Test case from http://se.mathworks.com/help/stats/quantile.html
            ([2, 4, 6, 8, 10, 12, 14], (4.5, 8.0, 11.5)),  # 7 = 4*1+3
            # Hand-made test case for 5 = 4*1+1 elements
            # n = 1 => q1 = 0.25*1st value + 0.75*2nd value
            #       => q3 = 0.75*4th value + 0.25*5th value
            ([3, 5, 7, 9, 11], (0.25*3 + 0.75*5, 7, 0.75*9+0.25*11)),
             ):
        got = compute_quartiles(test_data)
        if got != expected:
                raise AssertionError((got, expected))

aggregate_value_names = (
    "count",
    "avg", "std", "kurtosis", "skewness",
    "min", "q1", "median", "q3", "max",
    "paired_t",
    "p_value",
    )
    
def compute_aggregate_values(value_list):
    value_list = sorted(value_list)
    
    results = []

    # "count",
    n = len(value_list)
    assert n > 0
    results.append(n)
    
    # "avg", "std", "kurtosis", "skewness"
    avg = sum(value_list)/n
    results.append(avg)

    if n > 1:
        variance = online_variance(value_list)
        std = variance**0.5
    else:
        std = None
    results.append(std)

    if n > 2:
        kurtosis = online_kurtosis(value_list)
    else:
        kurtosis = None
    results.append(kurtosis)

    # This is 'sample skewness' from https://en.wikipedia.org/wiki/Skewness#Sample_skewness
    if n > 2:
        skew_top = sum((value-avg)**3 for value in value_list)/n
        skew_bot = (sum((value-avg)**2 for value in value_list)/(n-1))**1.5
        if skew_top:
            skewness = skew_top / skew_bot
        else:
            # 0/0 => 0.0
            skewness = 0.0
    else:
        skewness = None
    results.append(skewness)
    

    # "min", "q1", "median", "q3", "max",
    if n > 0:
        # XXX Should I use the average between the two middle points if there
        # is an even number of values?
        results.append(value_list[0])
        results.extend(compute_quartiles(value_list))
        results.append(value_list[-1])
    else:
        results.extend((None, None, None, None, None))

    
    # "paired_t",
    if n > 1:
        if std == 0.0:
            # MySQL doesn't handle infinity. Use 100000000 as the upper limit
            t = 100000000
        else:
            t = (avg / std) * n**0.5
            if t > 100000000:
                t = 100000000
    else:
        t = None
    results.append(t)

    # "p_value" from paired_t
    if n > 1:
        if std == 0.0:
            # XXX should I return this?
            p = None
        else:
            p = stats.t.sf(np.abs(t), n-1)*2
            # MySQL doesn't handle infinity. Use 100000000 as the upper limit
            if p > 100000000:
                p = 100000000
            
    else:
        p = None
    results.append(p)


    return results

def test_aggregate_values():
    values = compute_aggregate_values([2., 1., 3., 6.])
    assert len(aggregate_value_names) == len(values)
    count, avg, std, kurtosis, skewness, min_, q1, median, q3, max_, paired_t, p_value = values
    assert count == 4, count
    assert avg == 3.0, avg
    assert kurtosis == -1.0, kurtosis
    assert round(skewness, 3) == 0.446, skewness
    assert min_ == 1.0, min_
    assert q1 == 1.5, q1
    assert median == 2.5, median
    assert q3 == 4.5, q1
    assert max_ == 6.0, max_
    assert round(paired_t, 3) == 2.777, paired_t
    assert round(p_value, 3) == 0.069, p_value
test_aggregate_values()


class RuleEnvironment(object):
    __slots__ = ("idx", "property_value_lists")
    def __init__(self, idx, property_value_lists):
        self.idx = idx
        self.property_value_lists = property_value_lists

    def append_pair_properties(self, property_values1, property_values2):
        if not (property_values1 and property_values2):
            return
        if self.property_value_lists is None:
            self.property_value_lists = [[] for _ in property_values1]
            
        for property_list, value1, value2 in zip(self.property_value_lists,
                                                property_values1,
                                                property_values2):
            if value1 is None:
                continue
            if value2 is None:
                continue
            delta = value2-value1
            property_list.append(delta)
