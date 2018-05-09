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

import re

from rdkit import Chem
import itertools
from . import smiles_syntax # for validation

#####

class EnumerationLabel(object):
    NO_ENUMERATION = "N"
    CONSTANT_UP_ENUMERATION = "C"
    VARIABLE_UP_ENUMERATION = "V"

class Fragmentation(object):
    __slots__ = ("num_cuts", "enumeration_label",
                 "variable_num_heavies", "variable_symmetry_class", "variable_smiles",
                 "attachment_order",
                 "constant_num_heavies", "constant_symmetry_class", "constant_smiles", "constant_with_H_smiles")
    def __init__(self,
                 num_cuts, enumeration_label,
                 variable_num_heavies, variable_symmetry_class, variable_smiles,
                 attachment_order,
                 constant_num_heavies, constant_symmetry_class, constant_smiles, constant_with_H_smiles):
        self.num_cuts = num_cuts
        self.enumeration_label = enumeration_label
        self.variable_num_heavies = variable_num_heavies
        self.variable_symmetry_class = variable_symmetry_class
        self.variable_smiles = variable_smiles
        self.attachment_order = attachment_order
        self.constant_num_heavies = constant_num_heavies
        self.constant_symmetry_class = constant_symmetry_class
        self.constant_smiles = constant_smiles
        self.constant_with_H_smiles = constant_with_H_smiles

    def __repr__(self):
        return ("Fragmentation({self.num_cuts}, {self.enumeration_label!r}, "
                "{self.variable_num_heavies}, {self.variable_symmetry_class!r}, {self.variable_smiles!r}, "
                "{self.attachment_order!r}, "
                "{self.constant_num_heavies}, {self.constant_symmetry_class!r}, "
                "{self.constant_smiles!r}, {self.constant_with_H_smiles!r})").format(
                    self=self)

    def get_unique_key(self):
        return "%s.%s.%s" % (self.attachment_order, self.variable_smiles, self.constant_smiles)



#####

# TODO: Move some of these into smiles_syntax.py
    
# Extract just the atom terms (no closures)
_atom_pat = re.compile(r"""
(
 Cl? |
 Br? |
 [NOSPFIbcnosp] |
 \[[^]]*\]
)
""", re.X)
_atom_and_dot_disconnect_pat = re.compile(r"""
(
 Cl? |
 Br? |
 [NOSPFIbcnosp] |
 \[[^]]*\] |
 \* |
 \.
)
""", re.X)

def count_num_heavies(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def cansmiles(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def get_atom_order_in_smiles(mol):
    s = mol.GetProp("_smilesAtomOutputOrder")
    #print("Decode", s)
    positions = []
    assert s[0] == "[", s
    i = 1
    while 1:
        j = s.find(",", i)
        if j == -1:
            assert s[i:] == "]", (s, i, s[i:])
            break
        assert j > i, (s, i, s[i:])
        order = s[i:j]
        assert order.isdigit(), (s, i, j, order)
        i = j + 1
        positions.append(int(order))
    
    assert len(positions) == mol.GetNumAtoms(), (s, positions, mol.GetNumAtoms())
    assert len(set(positions)) == len(positions), positions
    return positions



def fragment_on_atom_pairs(mol, atom_pairs):
    bonds = []
    bond_dirs = {}
    dummy_labels = []
    
    label = 2
    for a1, a2 in atom_pairs:
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond.IsInRing():
            raise ValueError("Cannot fragment a ring bond (between %d and %d)"
                             % (a1, a2))
            
        bonds.append(bond.GetIdx())
        bond_dir = bond.GetBondDir()
        #print("bond is", bond_dir)

        if bond.GetBeginAtomIdx() == a1:
            dummy_labels.append( (label, label+1) )
            if bond_dir == Chem.BondDir.ENDDOWNRIGHT:
                bond_dirs[(a1, label+1)] = Chem.BondDir.ENDDOWNRIGHT
                bond_dirs[(a2, label  )] = Chem.BondDir.ENDUPRIGHT
            elif bond_dir == Chem.BondDir.ENDUPRIGHT:
                bond_dirs[(a1, label+1)] = Chem.BondDir.ENDUPRIGHT
                bond_dirs[(a2, label  )] = Chem.BondDir.ENDDOWNRIGHT
        else:
            # swapped
            dummy_labels.append( (label+1, label) )
            if bond_dir == Chem.BondDir.ENDUPRIGHT:
                bond_dirs[(a1, label+1)] = Chem.BondDir.ENDDOWNRIGHT
                bond_dirs[(a2, label  )] = Chem.BondDir.ENDUPRIGHT
            elif bond_dir == Chem.BondDir.ENDDOWNRIGHT:
                bond_dirs[(a1, label+1)] = Chem.BondDir.ENDUPRIGHT
                bond_dirs[(a2, label  )] = Chem.BondDir.ENDDOWNRIGHT
        
        label += 2

    ## print("bond_dirs", bond_dirs)
    ## print("bonds:", bonds)
    ## print("dummyLabels:", dummy_labels)
    fragmented_mol = Chem.FragmentOnBonds(mol, bonds, dummyLabels=dummy_labels)
    ## print("fragmented_mol:", cansmiles(fragmented_mol))

    dummy_pairs = [[] for _ in atom_pairs]
    for atom in fragmented_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            label = atom.GetIsotope()
            i = label // 2 - 1
            dummy_pairs[i].append(atom.GetIdx())
            atom.SetIsotope(0)

            for bond in atom.GetBonds():
                other_atom_id = bond.GetOtherAtomIdx(atom.GetIdx())
                #print("look for", (other_atom_id, label))
                bond_dir = bond_dirs.get( (other_atom_id, label), None)
                #print("set to", bond_dir)
                if bond_dir is not None:
                    bond.SetBondDir(bond_dir)
                break
            else:
                raise AssertionError
    
    other_atom_table = {}
    for a1, a2 in dummy_pairs:
        other_atom_table[a1] = a2
        other_atom_table[a2] = a1
    
    return fragmented_mol, other_atom_table #dummy_pairs

def get_num_heavies_from_smiles(smiles):
    num_atoms = 0
    for m in _atom_pat.finditer(smiles):
        text = m.group()
        if text == "[H]" or "*" in text:
            continue
        num_atoms += 1
    return num_atoms

def get_component_atom_symbols(smiles):
    components = []
    idx = 0
    component = []
    for m in _atom_and_dot_disconnect_pat.finditer(smiles):
        text = m.group()
        if text == ".":
            components.append(component)
            component = []
        else:
            component.append((idx, text))
            idx += 1
    components.append(component)
    return components
    

## def powerset(iterable):
##     "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
##     s = list(iterable)
##     return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


###

# I can't use Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
# because that does not identify ring stereocenters.
# The following is from Code/GraphMol/Chirality.cpp:isAtomPotentialChiralCenter.

def is_potential_chiral_center(atom, atom_ranks):
    d = atom.GetTotalDegree()
    if d not in (3, 4):
        # Only support tetrahedral chirality
        return False
    
    seen = set()
    atom_idx = atom.GetIdx()
    for bond in atom.GetBonds():
        other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
        rank = atom_ranks[other_atom_idx]
        if rank in seen:
            return False
        seen.add(rank)

    num_neighbors = len(seen)
    assert num_neighbors <= 4, num_neighbors
    if num_neighbors == 4:
        # 4 real neighbors (which means TotalNumHs is 0)
        # and all distinct.
        return True

    if num_neighbors < 3:
        # Not a stereocenter. Might even have 2 hydrogens.
        return False
    
    # 3 real neighbors. If there's a hydrogen, accept it.
    if atom.GetTotalNumHs() == 1:
        return True

    atomic_num = atom.GetAtomicNum()
    if atomic_num == 16 or atomic_num == 34:
        valence = atom.GetExplicitValence()
        if (valence == 4) or (valence == 3 and atom.GetFormalCharge() == 1):
            return True
    return False


CHIRAL_TAGS = (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
               Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

def get_chiral_flags(mol, atom_ranks):
    flags = []
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() in CHIRAL_TAGS:
            flg = 1
        elif is_potential_chiral_center(atom, atom_ranks):
            flg = 2
        else:
            flg = 0
        flags.append(flg)
    return flags

_H_cache = {}

def replace_wildcard_with_H(smiles):
    # The cache gives about 2% overall performance improvement.
    # My tests suggest there's about a 50% cache hit.
    try:
        return _H_cache[smiles]
    except KeyError:
        pass
    if smiles.count("[*]") == 1:
        smiles_with_H = smiles.replace("[*]", "[H]")
    elif smiles.count("*") == 1:
        smiles_with_H = smiles.replace("*", "[H]")
    else:
        raise AssertionError("Could not find the '*' atom")
        
    new_smiles = Chem.CanonSmiles(smiles_with_H)
    if len(_H_cache) > 10000:
        _H_cache.clear()
    _H_cache[smiles] = new_smiles
    return new_smiles

def make_single_cut(mol, atom_pair, chiral_flags):
    fragmented_mol, other_atom_table = fragment_on_atom_pairs(mol, [atom_pair])
    frag1_indices, frag2_indices = Chem.GetMolFrags(fragmented_mol)

    # Remove the indices for the wildcard atoms (should be the last two atoms in the molecule)
    num_atoms = fragmented_mol.GetNumAtoms()
    a1, a2 = num_atoms - 1, num_atoms - 2
    assert fragmented_mol.GetAtomWithIdx(a1).GetAtomicNum() == 0
    assert fragmented_mol.GetAtomWithIdx(a2).GetAtomicNum() == 0
    
    frag1_smiles = Chem.MolFragmentToSmiles(fragmented_mol, frag1_indices, isomericSmiles=True)
    frag2_smiles = Chem.MolFragmentToSmiles(fragmented_mol, frag2_indices, isomericSmiles=True)
    
    frag1_num_atoms = get_num_heavies_from_smiles(frag1_smiles)
    frag2_num_atoms = get_num_heavies_from_smiles(frag2_smiles)
    
    # Determine the symmetry of both parts
    fragmented_mol.UpdatePropertyCache(strict=False) # XXX magic; without it I get a RuntimeError

    # Need to clear chiral tags which are no longer relevant because the new
    # wildcards are symmetric. The canonical SMILES output is affected by an
    # atom's chiral tag, even if the output doesn't denote chirality for that
    # atom. I need to clear the tags to get a truly canonical output.
    # See https://sourceforge.net/p/rdkit/mailman/message/35420297/ , from Greg
    # Landrum, on 2016-10-11 05:39:12 titled "identify chiral atoms which
    # became achiral after fragmenting".
    Chem.AssignStereochemistry(fragmented_mol, cleanIt=True, force=True)
    
    #   "getNumImplicitHs() called without preceding call to calcImplicitValence()"
    new_atom_ranks = Chem.CanonicalRankAtoms(fragmented_mol, breakTies=False)
    ## print("new_atom_ranks:", list(new_atom_ranks))
    new_chiral_flags = get_chiral_flags(mol, new_atom_ranks)

    up_enumerations = []
    for frag_indices in (frag1_indices, frag2_indices):
        ## print("indices", frag_indices)
        ## print("chiral_flags", len(chiral_flags), chiral_flags)
        ## print("new_chiral_flags", len(new_chiral_flags), new_chiral_flags)

        frag_indices_without_wildcard = [a for a in frag1_indices if a < a2]
        chiral_indices = get_new_stereocenter_indices(
            frag_indices_without_wildcard, chiral_flags, new_chiral_flags)
        up_enumeration = set()
        for chiral_assignment in chiral_enumerate(chiral_indices):
            for (atom_index, chiral_tag) in chiral_assignment:
                fragmented_mol.GetAtomWithIdx(atom_index).SetChiralTag(chiral_tag)
            up_smiles = Chem.MolFragmentToSmiles(fragmented_mol, frag_indices, isomericSmiles=True)
            up_enumeration.add(up_smiles)
        up_enumerations.append(up_enumeration)

    frag1_up_enumerations, frag2_up_enumerations = up_enumerations
   
    # fragment 1 is the constant part and 2 is variable.
    for ((constant_num_atoms, constant_smiles, constant_up_enumerations,
          variable_num_atoms, variable_smiles, variable_up_enumerations)) in (
              (frag1_num_atoms, frag1_smiles, frag1_up_enumerations, frag2_num_atoms, frag2_smiles, frag2_up_enumerations),
              (frag2_num_atoms, frag2_smiles, frag2_up_enumerations, frag1_num_atoms, frag1_smiles, frag1_up_enumerations),              ):

        constant_smiles_with_H = replace_wildcard_with_H(constant_smiles)
        yield Fragmentation(1, EnumerationLabel.NO_ENUMERATION,
                            variable_num_atoms, "1", variable_smiles,
                            "0",
                            constant_num_atoms, "1", constant_smiles, constant_smiles_with_H)

        # up-enumeration in the constant part
        for constant_up_smiles in constant_up_enumerations:
            yield Fragmentation(1, EnumerationLabel.CONSTANT_UP_ENUMERATION,
                                variable_num_atoms, "1", variable_smiles,
                                "0",
                                constant_num_atoms, "1", constant_up_smiles, replace_wildcard_with_H(constant_up_smiles))
            
        # up-enumeration in the variable part
        for variable_up_smiles in variable_up_enumerations:
            yield Fragmentation(1, EnumerationLabel.VARIABLE_UP_ENUMERATION,
                                variable_num_atoms, "1", variable_up_smiles,
                                "0",
                                constant_num_atoms, "1", constant_up_smiles, constant_smiles_with_H)
    


def _get_bonds_from_atom_pairs(mol, atom_pairs):
    bonds = []
    for a1, a2 in atom_pairs:
        bonds.append(mol.GetBondBetweenAtoms(a1, a2))
    return bonds

def _get_variable_index(smiles_list):
    num_cuts = len(smiles_list) - 1
    for i, smiles in enumerate(smiles_list):
        n = smiles.count("*")
        if n == 1:
            continue
        if n != num_cuts:
            # 3 cuts but not on a central core
            assert n == 2 and num_cuts == 3, (smiles_list, i)
            return None
        return i
    raise AssertionError(smiles_list)

def get_symmetry_class(a, b, c=None):
    if c is None:
        if a == b:
            return "11"
        else:
            return "12"
    if a == b:
        if b == c:
            return "111"
        return "112"
    if a == c:
        return "121"
    if b == c:
        return "122"
    return "123"
    

def _init_canonical_order():
    canonical_order = {}

    def get_connection(symm_group1, symm_group2, perm):
        terms = []
        for i, p in enumerate(perm):
            j = int(p)
            terms.append(symm_group1[i] + symm_group2[j])
        terms.sort()
        return terms
        
    # Many permutations can give the same mapping between
    # constant and variable location. The canonical mapping
    # is the lexically smallest of the possible mappings.
    for symmetry_groups, permutations in (
            (("11", "12"), ("01", "10")),
            (("111", "112", "121", "122", "123"),
             ("012", "021", "102", "120", "201", "210")),
             ):
        ordered_permutations = sorted(permutations)
        
        for symm_group1 in symmetry_groups:
            for symm_group2 in symmetry_groups:
                for perm in permutations:
                    target_connection = get_connection(
                        symm_group1, symm_group2, perm)
                    for canonical_perm in ordered_permutations:
                        if get_connection(
                            symm_group1, symm_group2, canonical_perm) == target_connection:
                            canonical_order[symm_group1, symm_group2, perm] = canonical_perm
                            break
                    else:
                        raise AssertionError
    return canonical_order

CANONICAL_ATTACHMENT_ORDER = _init_canonical_order()

def get_chiral_difference(atom_indices, old_chiral_flags, new_chiral_flags):
    num_chirals = num_lost_chirals = num_new_stereocenters = 0
    for atom_index in atom_indices:
        old_flag = old_chiral_flags[atom_index]
        new_flag = new_chiral_flags[atom_index]
        if old_flag == 0:
            if new_flag == 0:
                pass # expected
            elif new_flag == 1:
                raise AssertionError("that shouldn't happen")
            else:
                num_new_stereocenters += 1
        elif old_flag == 1:
            if new_flag == 0:
                raise AssertionError("that also shouldn't happen")
            elif new_flag == 1:
                num_chirals += 1 # still chiral
            else:
                raise AssertionError("That was unexpected")
        elif old_flag == 2:
            if new_flag == 0:
                # Wasn't chiral, could have been, but now cannot.
                pass
            elif new_flag == 1:
                raise AssertionError("chiral was *added*?")
    return num_chirals, num_lost_chirals, num_new_stereocenters

# Atoms which weren't a stereocenter due to symmetry but which, after
# fragmentation, can be a stereocenter
def up_enumerate(fragmented_mol, constant_atom_indices, variable_atom_indices,
                 chiral_flags, new_chiral_flags):
    yield EnumerationLabel.NO_ENUMERATION, None
    constant_indices = get_new_stereocenter_indices(constant_atom_indices, chiral_flags, new_chiral_flags)
    #print("test", constant_atom_indices)
    #print([(i, chiral_flags[i], new_chiral_flags[i]) for i in constant_atom_indices])
    #print("constant_indices", constant_indices)
    
    if constant_indices:
        for chiral_enumeration in chiral_enumerate(constant_indices):
            yield EnumerationLabel.CONSTANT_UP_ENUMERATION, chiral_enumeration
    
    variable_indices = get_new_stereocenter_indices(variable_atom_indices, chiral_flags, new_chiral_flags)
    if variable_indices:
        for chiral_enumeration in chiral_enumerate(variable_indices):
            yield EnumerationLabel.VARIABLE_UP_ENUMERATION, chiral_enumeration

def get_new_stereocenter_indices(atom_indices, old_chiral_flags, new_chiral_flags):
    stereocenter_indices = []
    for atom_index in atom_indices:
        old_flag = old_chiral_flags[atom_index]
        new_flag = new_chiral_flags[atom_index]
        if old_flag == 0 and new_flag == 2:
            stereocenter_indices.append(atom_index)
    return stereocenter_indices

def chiral_enumerate(indices):
    chiral_tags = (Chem.CHI_UNSPECIFIED,
                   Chem.CHI_TETRAHEDRAL_CW,
                   Chem.CHI_TETRAHEDRAL_CCW)
    terms = []
    for index in indices:
        terms.append((index, chiral_tag) for chiral_tag in chiral_tags)
    it = itertools.product(*terms)
    next(it)  # The first one is the input structure
    return it

        
def make_multiple_cuts(mol, atom_pairs, chiral_flags):
    num_cuts = len(atom_pairs)
    assert num_cuts >= 2, num_cuts
    fragmented_mol, other_atom_table = fragment_on_atom_pairs(mol, atom_pairs)

    # Figure out which atoms are in the variable part and which atoms are in the constant part.

    constant_atom_indices = []
    variable_atom_indices = []
    for atom_indices in Chem.GetMolFrags(fragmented_mol):
        non_wildcard_indices = []
        for atom_index in atom_indices:
            if fragmented_mol.GetAtomWithIdx(atom_index).GetAtomicNum() != 0:
                non_wildcard_indices.append(atom_index)
        num_wildcard_atoms = len(atom_indices) - len(non_wildcard_indices)
        if num_wildcard_atoms == 1:
            constant_atom_indices.extend(non_wildcard_indices)
        elif num_wildcard_atoms == num_cuts:
            variable_atom_indices.extend(non_wildcard_indices)
        else:
            # Did not cut into core+rgroups
            return

    # Determine the symmetry of the variable part
    fragmented_mol.UpdatePropertyCache(strict=False) # XXX magic; without it I get a RuntimeError
    Chem.AssignStereochemistry(fragmented_mol, cleanIt=True, force=True)

    #   "getNumImplicitHs() called without preceding call to calcImplicitValence()"
    new_atom_ranks = Chem.CanonicalRankAtoms(fragmented_mol, breakTies=False)
    new_chiral_flags = get_chiral_flags(mol, new_atom_ranks)

    seen_smiles = set()
    #     
    for enumeration_label, chiral_assignments in up_enumerate(
            fragmented_mol, constant_atom_indices, variable_atom_indices, chiral_flags, new_chiral_flags):
        if enumeration_label == EnumerationLabel.NO_ENUMERATION:
            assert chiral_assignments is None
            atom_ranks = new_atom_ranks
            ## print("reused:", list(atom_ranks))
        else:
            for (atom_index, chiral_tag) in chiral_assignments:
                fragmented_mol.GetAtomWithIdx(atom_index).SetChiralTag(chiral_tag)
            fragmented_mol.ClearComputedProps() # XXX Do I need this?
            atom_ranks = Chem.CanonicalRankAtoms(fragmented_mol, breakTies=False)
            ## print("computed:", list(atom_ranks))


        # Work in SMILES space so we find a canonical mapping between the
        # unlabeled canonical variable and canonical constant parts.
        smiles = cansmiles(fragmented_mol)
        #print("smiles", smiles)
        
        # The up-enumeration may have several ways to generate the same structure.
        # For example, flipping two "@"s to "@@"s may leave the structure unchanged.
        if smiles in seen_smiles:
            continue
        seen_smiles.add(smiles)
    
        # Figure out which is the variable/core structure.
        # It's the one with the most "*"s on it (must equal the number of cuts)
        frag_smiles_list = smiles.split(".")
        assert len(frag_smiles_list) == num_cuts+1, smiles
        variable_component_index = _get_variable_index(frag_smiles_list)
        if variable_component_index is None:
            # 3 cuts but no fragment with three "*"s
            raise AssertionError(("I already checked for this", smiles))
    
        #print("core is at", variable_component_index)

        # Get the mapping from position in the SMILES string to atom index in the molecule
        smiles_index_to_atom_index = get_atom_order_in_smiles(fragmented_mol)

        # Determine the constant part (the rgroups)
        constant_component_indices = list(range(num_cuts+1))
        del constant_component_indices[variable_component_index]
        constant_smiles_list = [frag_smiles_list[i] for i in constant_component_indices]
        assert len(constant_smiles_list) == num_cuts
        
        # Find the connection points on the variable part
        component_atom_symbols = get_component_atom_symbols(smiles)
        variable_connection_atom_indices = []
        variable_atom_indices2 = []
        for smiles_index, smiles_symbol in component_atom_symbols[variable_component_index]:
            atom_index = smiles_index_to_atom_index[smiles_index]
            if "*" in smiles_symbol:
                variable_connection_atom_indices.append(atom_index)
            else:
                variable_atom_indices2.append(atom_index) # XXX Remove
        assert sorted(variable_atom_indices) == sorted(variable_atom_indices2), (
            sorted(variable_atom_indices), sorted(variable_atom_indices2))
            
        assert len(variable_connection_atom_indices) == num_cuts
        
        #print("variable_connection_atom_indices", variable_connection_atom_indices)
        variable_symmetry_class = get_symmetry_class(*(
            atom_ranks[atom_index] for atom_index in variable_connection_atom_indices))
        
        # Determine the symmetry of the constant part (the rgroups)

        constant_symmetry_class = get_symmetry_class(*constant_smiles_list)

        # Figure out which R-groups in the constant part correspond to the
        # attachment points in the core/variable part.
        atom_index_to_rgroup_label = {}
        constant_atom_indices = []
        for rgroup_id, component_i in enumerate(constant_component_indices):
            rgroup_label = str(rgroup_id)
            for (smiles_index, smiles_symbol) in component_atom_symbols[component_i]:
                atom_index = smiles_index_to_atom_index[smiles_index]
                atom_index_to_rgroup_label[atom_index] = rgroup_label
                if "*" not in smiles_symbol:
                    constant_atom_indices.append(atom_index)

        attachment_order = "".join(atom_index_to_rgroup_label[other_atom_table[atom_index]]
                                     for atom_index in variable_connection_atom_indices)
        # Figure the canonical attachment order
        canonical_attachment_order = CANONICAL_ATTACHMENT_ORDER[
            variable_symmetry_class,
            constant_symmetry_class,
            attachment_order]

        # Figure out which atoms in the variable part are still chiral
        ## fragmented_chiral_flags = get_chiral_flags(fragmented_mol, atom_ranks)
        ## variable_num_chirals, variable_num_lost_chirals, variable_num_new_stereocenters = \
        ##   get_chiral_difference(variable_atom_indices2, chiral_flags, fragmented_chiral_flags)

        ## constant_num_chirals, constant_num_lost_chirals, constant_num_new_stereocenters = \
        ##   get_chiral_difference(constant_atom_indices2, chiral_flags, fragmented_chiral_flags)

        variable_smiles = frag_smiles_list[variable_component_index]
        constant_smiles = ".".join(constant_smiles_list)
        ## print("variable_smiles:", variable_smiles)
        ## print("constant_smiles:", constant_smiles)

        # Test that I can reconnect
        if 0:
            offsets = [int(c) for c in canonical_attachment_order]
            var_part = smiles_syntax.convert_wildcards_to_closures(variable_smiles, offsets)
            const_part = smiles_syntax.convert_wildcards_to_closures(constant_smiles, list(range(num_cuts)))
            smi = Chem.CanonSmiles(var_part + "." + const_part, 0)
            expected_smi = Chem.MolToSmiles(mol)
            if smi != expected_smi:
                print("     Got:", smi)
                print("Expected:", expected_smi)
            assert smi == expected_smi, (smi, expected_smi)

        ## print("Fragmentation")
        ## print(get_num_heavies_from_smiles(variable_smiles), variable_symmetry_class, variable_smiles)

        yield Fragmentation(
            num_cuts,
            enumeration_label,
            get_num_heavies_from_smiles(variable_smiles), variable_symmetry_class, variable_smiles,
            canonical_attachment_order,
            get_num_heavies_from_smiles(constant_smiles), constant_symmetry_class, constant_smiles, None,
            )
    
            
def fragment_mol(mol, fragment_filter, num_heavies=None):
    try:
        for x in _fragment_mol(mol, fragment_filter, num_heavies):
            yield x
    except:
        ## import traceback
        ## traceback.print_exc()
        raise
    
def _fragment_mol(mol, fragment_filter, num_heavies=None):
    
    cut_lists = fragment_filter.get_cut_lists(mol)

    if not cut_lists:
        return

    seen = set()
    
    if num_heavies is None:
        num_heavies = count_num_heavies(mol)

    # Identify atoms that are chiral (assigned and unassigned)in parent compound
    # 0 means not chiral, 1 means assigned, 2 means unassigned
    atom_ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    chiral_flags = get_chiral_flags(mol, atom_ranks)

    for cut_list in cut_lists:
        num_cuts = len(cut_list)
        #print("num_cuts", num_cuts)
        if num_cuts == 1:
            fragmentations = make_single_cut(mol, cut_list[0], chiral_flags)
        else:
            fragmentations = make_multiple_cuts(mol, cut_list, chiral_flags)
        
        for fragmentation in fragmentations:
            key = fragmentation.get_unique_key() # XXX + "012" + YYY
            if key not in seen:
                seen.add(key)
                yield fragmentation


### fragment on hydrogens

# NOTE: this is hard-coded to match the string used in
# index_algorithm.py's load_fragment_index
_hydrogen_cut_smiles = "[*][H]"


# (Used by the 'transform' code.) Given a molecule, synthesize fragments
# based on breaking a hydrogen.
def get_hydrogen_fragmentations(smiles, num_heavies):
    fragmentations = []
    mol = Chem.MolFromSmiles(smiles)
    seen = set()
    for atom in mol.GetAtoms():
        # All hydrogens are equivalent, so only need 1 h-fragment per atoms
        if atom.GetNumImplicitHs() > 0:
            emol = Chem.EditableMol(mol)
            # Add the "*", single-bonded to the atom
            wildcard_atom_idx = emol.AddAtom(Chem.Atom(0))
            emol.AddBond(atom.GetIdx(), wildcard_atom_idx, Chem.BondType.SINGLE)
            cut_mol = emol.GetMol()
            cut_smiles = Chem.MolToSmiles(cut_mol, isomericSmiles=True)
            if cut_smiles in seen:
                continue
            seen.add(cut_smiles)
            new_fragmentation = Fragmentation(
                1, EnumerationLabel.NO_ENUMERATION,
                0, "1", _hydrogen_cut_smiles,
                "0",
                num_heavies, "1", cut_smiles, smiles)
            fragmentations.append(new_fragmentation)
    return fragmentations


###

# If there is an explicit '[H]' in the SMILES then fragment only on that.

_hydrogen_cut_pat = Chem.MolFromSmarts("[!#1]-[0#1v1!+!-]")

def fragment_molecule_on_explicit_hydrogens(smiles):
    num_heavies = get_num_heavies_from_smiles(smiles)
    smiles_with_H = Chem.CanonSmiles(smiles)
    input_mol = Chem.MolFromSmiles(smiles, sanitize=False) # use santize=False to preserve explicit hydrogens
    Chem.SanitizeMol(input_mol, Chem.SANITIZE_ALL)
    
    cut_pairs = input_mol.GetSubstructMatches(_hydrogen_cut_pat)

    fragmentations = []
    for cut_pair in cut_pairs:
        bond_idx = input_mol.GetBondBetweenAtoms(*cut_pair).GetIdx()
        fragmented_mol = Chem.FragmentOnBonds(input_mol, [bond_idx], dummyLabels=[(0, 0)])
        new_smiles = Chem.MolToSmiles(fragmented_mol, isomericSmiles=True)
        
        left, mid, right = new_smiles.partition(".")
        assert mid == ".", new_smiles

        if left == "[*][H]":  # Hard-coded
            cut_smiles = right
        elif right == "[*][H]":
            cut_smiles = left
        else:
            raise AssertionError("did not split hydrogen correctly: %r %r"
                                  % (smiles, new_smiles))

        if "[H]" in cut_smiles:
            # If there were multiple [H] atoms, then we cut on one but others remain.
            # Recanonicalize to remove them.
            cut_smiles = Chem.CanonSmiles(cut_smiles)

        new_fragmentation = Fragmentation(
            1, EnumerationLabel.NO_ENUMERATION,
            0, "1", "[*][H]",
            "0",
            num_heavies, "1", cut_smiles,
            None)

        fragmentations.append(new_fragmentation)

    return fragmentations
        

#, dummyLabels=[(1, 1)]
