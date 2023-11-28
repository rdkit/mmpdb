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

# NOTE: There is configuration information in three files!
# 1) fragment_types.py -- the data types in this file
# 2) fragment_schema.sql -- defines the SQL schema
# 3) fragment_db.py -- defines the mapping from SQL to the data types

# I tried using SQLAlchemy to merge them into one but that
# added a ~30% overhead to fragment generation.

import dataclasses
from dataclasses import dataclass
from typing import List, Optional

try:
    from typing import Literal
except ImportError:
    # Compatability for Python 3.7
    # ("typing.Literal" was added in Python 3.8)
    # Make Literai[x] return Any
    from typing import Any

    class LiteralClass(dict):
        def __missing__(self, x):
            return Any

    Literal = LiteralClass()

class FragmentationFailure(Exception):
    pass

@dataclass
class FragmentOptions:
    max_heavies: Optional[int]
    max_rotatable_bonds: Optional[int]
    rotatable_smarts: str
    cut_smarts: str

    num_cuts: Literal[1, 2, 3]
    method: Literal["chiral"]
    salt_remover: str
    min_heavies_per_const_frag: int
    max_up_enumerations: int
    min_heavies_total_const_frag: int = 0 # added in fragment format v4, 28 Nov 2023

    def to_dict(self):
        return dataclasses.asdict(self)

    ## # For backwards compatability (still needed?) -- TODO: remove
    ## def to_text_settings(self):
    ##     def _none(x):
    ##         return "none" if x is None else str(x)
    ##     return tuple((name, _none(value))
    ##                      for (name, value) in dataclasses.asdict(self).items())

    def get_fragment_filter(self):
        return get_fragment_filter(self)


@dataclass
class Fragmentation:
    num_cuts: int
    enumeration_label: str
    variable_num_heavies: int
    variable_symmetry_class: str
    variable_smiles: str
    attachment_order: str
    constant_num_heavies: int
    constant_symmetry_class: str
    constant_smiles: str
    constant_with_H_smiles: Optional[str]

    def get_unique_key(self):
        return "%s.%s.%s" % (
            self.attachment_order,
            self.variable_smiles,
            self.constant_smiles,
        )


@dataclass
class FragmentRecord:
    """An input structure record which could be parsed and fragmented, and its fragmentations"""

    id: str
    input_smiles: str
    num_normalized_heavies: int
    normalized_smiles: str
    fragmentations: List[Fragmentation]
    errmsg = None  # make it easier to tell if this is an error record


@dataclass
class FragmentErrorRecord(object):
    """An input structure record which could not be parsed and fragmented"""

    id: str
    input_smiles: str
    errmsg: str


## Exceptions


class FragmentValueError(ValueError):
    def __init__(self, name, value, reason):
        self.name = name
        self.value = value
        self.reason = reason

    def __str__(self):
        return "Error with %s (%r): %s" % (self.name, self.value, self.reason)

    def __repr__(self):
        return "FragmentValueError(%r, %r, %r)" % (self.name, self.value, self.reason)


class FragmentFormatError(ValueError):
    def __init__(self, reason, location):
        self.reason = reason
        self.location = location

    def __repr__(self):
        return "FragmentFormatError(%r, %r)" % (self.reason, self.location)

    def __str__(self):
        return "%s, %s" % (self.reason, self.location.where())


## FragmentFilter


def parse_rotatable_smarts(smarts):
    from rdkit import Chem

    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError("unable to parse SMARTS")
    if pattern.GetNumAtoms() != 2:
        raise ValueError("rotatable SMARTS must match exactly two atoms")
    if pattern.GetNumBonds() != 1:
        raise ValueError("rotatable SMARTS must connect both atoms")

    return pattern


def parse_cut_smarts(smarts):
    from rdkit import Chem
    from . import smarts_aliases

    if smarts in smarts_aliases.cut_smarts_aliases_by_name:
        smarts = smarts_aliases.cut_smarts_aliases_by_name[smarts].smarts
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError("unable to parse SMARTS")
    if pattern.GetNumAtoms() != 2:
        raise ValueError("cut SMARTS must match exactly two atoms")
    if pattern.GetNumBonds() != 1:
        raise ValueError("cut SMARTS must connect both atoms")
    return pattern


def parse_salt_remover(salt_remover_filename):
    from rdkit.Chem import SaltRemover

    if salt_remover_filename == "<none>":
        # Don't use a salt remover
        return None
    elif salt_remover_filename == "<default>":
        remover = SaltRemover.SaltRemover()
    else:
        # raises an I/O error on failure
        remover = SaltRemover.SaltRemover(salt_remover_filename)

    # Work-around rdkit issue #541, where SMARTS failures get
    # turned into None patterns instead of raising an exception
    # during SaltRemover initialization.
    for i, pattern in enumerate(remover.salts, 1):
        if pattern is None:
            raise ValueError("SMARTS pattern #%d is invalid" % (i,))
    return remover


def parse_method(method):
    from . import fragment_algorithm

    if method in ("", None, "chiral"):
        return fragment_algorithm.fragment_mol
    raise ValueError("must be 'chiral'")


def parse_num_cuts(num_cuts):
    if num_cuts in (1, 2, 3):
        return num_cuts
    raise ValueError("must be 1, 2, or 3")


class FragmentFilter(object):
    def __init__(
        self,
        *,
        max_heavies,
        max_rotatable_bonds,
        rotatable_pattern,
        salt_remover,
        cut_pattern,
        num_cuts,
        method,
        options,
        min_heavies_per_const_frag,
        min_heavies_total_const_frag,
        max_up_enumerations,
    ):
        if rotatable_pattern is None:
            max_rotatable_bonds = None

        self.max_heavies = max_heavies
        self.max_rotatable_bonds = max_rotatable_bonds
        self.rotatable_pattern = rotatable_pattern
        self.salt_remover = salt_remover
        self.cut_pattern = cut_pattern
        self.num_cuts = num_cuts
        self.method = method
        self.options = options
        self.min_heavies_per_const_frag = min_heavies_per_const_frag
        self.min_heavies_total_const_frag = min_heavies_total_const_frag
        self.max_up_enumerations = max_up_enumerations

    def normalize(self, mol):
        # XXX Remove existing isotope labels?
        if self.salt_remover is not None:
            desalted_mol = self.salt_remover.StripMol(mol)
            if not mol.GetNumAtoms():
                return ("no non-salts", desalted_mol)
        else:
            desalted_mol = mol

        return None, desalted_mol

    def apply_filters(self, mol):
        num_heavies = 0
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num <= 1:
                if atomic_num == 0:
                    return "contains '*' atoms"
                else:
                    # should I worry about explicit hydrogens?
                    pass
            else:
                num_heavies += 1

        if num_heavies < 2:
            return "not enough heavy atoms"

        max_heavies = self.max_heavies
        if max_heavies is not None and num_heavies > max_heavies:
            return "too many heavy atoms"

        max_rotatable_bonds = self.max_rotatable_bonds
        if max_rotatable_bonds is not None:
            # The maxMatches (in RDKit at least up to early 2016) specifies
            # the number of non-unique matches. After those N are found, they
            # are filtered to remove duplicates. It's possible that up to
            # 1/2 of the bonds will be removed, so the limit must be twice
            # the max requested, plus 1
            matches = mol.GetSubstructMatches(
                self.rotatable_pattern,
                uniquify=True,
                maxMatches=max_rotatable_bonds * 2 + 1,
            )
            if len(matches) > max_rotatable_bonds:
                return "too many rotatable bonds"

        return None

    def get_cut_atom_pairs(self, mol):
        return mol.GetSubstructMatches(self.cut_pattern, uniquify=True)

    def get_cut_lists(self, mol):
        atom_pairs = self.get_cut_atom_pairs(mol)

        # Generate all of the atoms lists with 1, 2, or 3 cuts
        cut_lists = []
        for i, first_pair in enumerate(atom_pairs):
            cut_lists.append([first_pair])

            if self.num_cuts >= 2:
                for j, second_pair in enumerate(atom_pairs[i + 1 :], i + 1):
                    cut_lists.append([first_pair, second_pair])

                    if self.num_cuts >= 3:
                        for k, third_pair in enumerate(atom_pairs[j + 1 :], j + 1):
                            cut_lists.append([first_pair, second_pair, third_pair])
        return cut_lists


def get_fragment_filter(fragment_options):
    options = fragment_options

    def call(parse, name):
        value = getattr(options, name)
        try:
            return parse(value)
        except ValueError as err:
            raise FragmentValueError(name, value, str(err))

    max_heavies = options.max_heavies
    max_rotatable_bonds = options.max_rotatable_bonds
    rotatable_pattern = call(parse_rotatable_smarts, "rotatable_smarts")
    cut_pattern = call(parse_cut_smarts, "cut_smarts")

    num_cuts = call(parse_num_cuts, "num_cuts")
    method = call(parse_method, "method")
    min_heavies_per_const_frag = options.min_heavies_per_const_frag
    min_heavies_total_const_frag = options.min_heavies_total_const_frag
    try:
        salt_remover = call(parse_salt_remover, "salt_remover")
    except IOError as err:
        raise FragmentValueError("salt_remover", options.salt_remover, f"Cannot open salt file: {err}")
    
    max_up_enumerations = options.max_up_enumerations
    
    return FragmentFilter(
        max_heavies=max_heavies,
        max_rotatable_bonds=max_rotatable_bonds,
        rotatable_pattern=rotatable_pattern,
        salt_remover=salt_remover,
        cut_pattern=cut_pattern,
        num_cuts=num_cuts,
        method=method,
        options=fragment_options,
        min_heavies_per_const_frag=min_heavies_per_const_frag,
        min_heavies_total_const_frag=min_heavies_total_const_frag,
        max_up_enumerations=max_up_enumerations,
    )


####
