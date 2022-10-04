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

# Determine the environment around the attachment point(s) of the constant
# part using circular fingerprints.

import sys
import hashlib
import binascii

from rdkit import Chem
from rdkit.Chem import AllChem

if sys.version_info.major == 2:

    def fingerprint_to_text(s):
        return s


else:

    def fingerprint_to_text(s):
        return s.decode("ascii")


# Internal helper class
class EnvironmentCenters(object):
    def __init__(self, mol, atom_ids):
        self.mol = mol
        self.atom_ids = atom_ids

    def __repr__(self):
        return "EnvironmentCenters(%r, %r)" % (self.mol, self.atom_ids)


class EnvironmentFingerprint(object):
    __slots__ = ("radius", "smarts", "pseudosmiles", "parent_smarts")

    def __init__(self, radius, smarts, pseudosmiles, parent_smarts):
        self.radius = radius
        self.smarts = smarts
        self.pseudosmiles = pseudosmiles
        self.parent_smarts = parent_smarts
        # What about parent?

    def __repr__(self):
        return (
            f"EnvironmentFingerprint({self.radius}, {self.smarts!r}, "
            f"{self.pseudosmiles!r}, {self.parent_smarts})"
            )


def find_centers(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Could not parse the context SMILES %r" % (smiles,))

    # Figure out which atoms to use as the centers
    centers = {}
    for atom_idx, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() != 0:
            continue
        mapno = atom.GetIntProp("molAtomMapNumber")

        # Need to clear this because it can affect MolToSmiles() order.
        atom.ClearProp("molAtomMapNumber")
        assert mapno not in centers, (mapno, smiles)

        centers[mapno] = atom_idx

    if not centers:
        raise ValueError("No attachment points found in context SMILES %r" % (smiles,))

    n = len(centers)
    if n > 3:
        raise ValueError("Too many attachment points found in context SMILES %r" % (smiles,))
    if n == 3:
        try:
            ordered_centers = centers[1], centers[2], centers[3]
        except KeyError:
            raise ValueError("Expecting attachment points *:1, *:2, and *:3 in context SMILES %r" % (smiles,))
    elif n == 2:
        try:
            ordered_centers = centers[1], centers[2]
        except KeyError:
            raise ValueError("Expecting attachment points *:1 and *:2 in context SMILES %r" % (smiles,))
    elif n == 1:
        try:
            ordered_centers = centers[1],  # this is a one element tuple
        except KeyError:
            raise ValueError("Expecting attachment point *:1 in context SMILES %r" % (smiles,))

    return EnvironmentCenters(mol, ordered_centers)


# Get each of the atoms counts radii as it's being computed.
def iter_num_atoms_for_radii(centers, min_radius, max_radius):
    return _iter_num_atoms_for_radii(centers.mol, min_radius, max_radius, centers.atom_ids)


# My thought is to use this for testing.
def _iter_num_atoms_for_radii(mol, min_radius, max_radius, start_atoms):
    unique_atoms = set(start_atoms)
    assert len(start_atoms) == len(unique_atoms), "duplicate start atom"
    ignore_atoms = set(a for a in start_atoms if not is_heavy_atom(mol.GetAtomWithIdx(a)))

    yield (len(unique_atoms) - len(ignore_atoms))

    border_atoms = unique_atoms.copy()

    for radius in range(min_radius, max_radius):  # up to and including max_radius
        new_atoms = set()

        for atom in border_atoms:
            for neighbor_atom_obj in mol.GetAtomWithIdx(atom).GetNeighbors():
                neighbor_atom = neighbor_atom_obj.GetIdx()
                if neighbor_atom not in unique_atoms:
                    unique_atoms.add(neighbor_atom)
                    new_atoms.add(neighbor_atom)

                    if not is_heavy_atom(neighbor_atom_obj):
                        ignore_atoms.add(neighbor_atom)

        border_atoms = new_atoms
        yield (len(unique_atoms) - len(ignore_atoms))


def is_heavy_atom(atom_obj):
    # Exclude '*' and '[H]' atoms (but count '[1H]', '[2H]', etc.)
    eleno = atom_obj.GetAtomicNum()
    if eleno == 0 or (eleno == 1 and atom_obj.GetIsotope() == 0):
        return False
    return True



#### The Morgan atom connectivity invariants are:
#
## components.push_back(atom->getAtomicNum());
## components.push_back(atom->getTotalDegree());
## components.push_back(atom->getTotalNumHs());
## components.push_back(atom->getFormalCharge());
## int deltaMass = static_cast<int>(
##     atom->getMass() -
##     PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
## components.push_back(deltaMass);
#
# We need a SMARTS representation which matches that.
#  [#7,X3,H2,+1,R]
# This could be post-processed into a SMILES.

_aromatic_elements = (5, 6, 7, 8, 15, 16, 33, 34, 52)

def get_atom_symbols(mol, center_ids):
    aromatic_elements = _aromatic_elements
    
    atom_symbols = []
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        is_in_ring = atom.IsInRing()
        if atomic_num == 0:
            # '*' in SMARTS matches any atom.
            # Need to match only the atom with element number 0.
            element_symbol = "#0"
        elif (is_in_ring and atomic_num in aromatic_elements):
            # Unable to decide
            element_symbol = "#%d" % atomic_num
        else:
            # Must be aliphatic
            element_symbol = atom.GetSymbol()
            
        atom_symbols.append("[%s;X%d;H%d;%+d;%s]" % (
            element_symbol,
            atom.GetTotalDegree(),
            atom.GetTotalNumHs(),
            atom.GetFormalCharge(),
            ("!R", "R")[is_in_ring]))
            
    for attachment_point, center_id in enumerate(center_ids, 1):
        atom = mol.GetAtomWithIdx(center_id)
        atom_symbols[center_id] = "%s:%s]" % (
            atom_symbols[center_id][:-1], attachment_point
            )
    return atom_symbols

#### The Morgan bond connectivity invariants are:
#
##      bondInvariant = static_cast<int32_t>(bond->getBondType());
#
_bond_smarts_symbols = {
    Chem.BondType.SINGLE: "-",
    Chem.BondType.DOUBLE: "=",
    Chem.BondType.TRIPLE: "#",
    Chem.BondType.AROMATIC: ":",
    }

def get_bond_symbols(mol):
    return [_bond_smarts_symbols[bond.GetBondType()] for bond in mol.GetBonds()]

def get_environment_smarts_list_for_center(centers, max_radius, min_radius=0):
    assert max_radius >= min_radius, (min_radius, max_radius)
    assert min_radius >= 0, min_radius
    
    mol = centers.mol
    center_ids = centers.atom_ids
    
    atom_symbols = get_atom_symbols(mol, center_ids)
    bond_symbols = get_bond_symbols(mol)
    
    smarts_list = []

    atom_ids = set(center_ids)
    bond_ids = set()
    edge_atom_ids = atom_ids.copy()

    if min_radius <= 0:
        smarts = Chem.MolFragmentToSmiles(mol, list(atom_ids),
                                              atomSymbols = atom_symbols,
                                              bondSymbols = bond_symbols,
                                              isomericSmiles = False)
        smarts_list.append(smarts)

    # The Morgan fingerprint don't look at atom aromaticity.
    # Make everything aliphatic so MolFragmentToSmiles() doesn't
    # include atom aromaticity in its output.
    aromatics = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            atom.SetIsAromatic(False)
            aromatics.append(atom)
    
    for r in range(1, max_radius+1):
        new_edge_atom_ids = []
        for atom_id in edge_atom_ids:
            atom = mol.GetAtomWithIdx(atom_id)
            for bond in atom.GetBonds():
                bond_id = bond.GetIdx()
                if bond_id in bond_ids:
                    continue
                bond_ids.add(bond_id)

                other_atom_id = bond.GetOtherAtomIdx(atom_id)
                if other_atom_id in atom_ids:
                    continue
                atom_ids.add(other_atom_id)
                new_edge_atom_ids.append(other_atom_id)

        if r >= min_radius:
            smarts = Chem.MolFragmentToSmiles(mol, list(atom_ids), list(bond_ids),
                                                  atomSymbols = atom_symbols,
                                                  bondSymbols = bond_symbols,
                                                  isomericSmiles = False)

            # MolFragmentToSmiles() does not know the connection order,
            # even though the atom symbol includes the :1, :2, and :3.
            # Need to sort to put things in order.
            def get_term_key(term):
                if ":3]" in term: return 3
                if ":2]" in term: return 2
                if ":1]" in term: return 1
                raise AssertionError((term, smarts))
            terms = smarts.split(".")
            terms.sort(key = get_term_key)
            smarts = ".".join(terms)

            smarts_list.append(smarts)
        
        edge_atom_ids = new_edge_atom_ids

    # Revert aromaticity
    for atom in aromatics:
        atom.SetIsAromatic(True)
        
    return smarts_list

# 
# [#6;X3;H0;+0;R]-[C;X3;H0;+0;!R](=[C;X3;H1;+0;!R])-[*;X1;H0;+0;!R:2].[C;X4;H3;+0;!R]-[*;X1;H0;+0;!R:1]
#

import re
pt = Chem.GetPeriodicTable()
_aliphatic_symbol_lookup = dict((("#%d" % i), pt.GetElementSymbol(i)) for i in _aromatic_elements)
_aromatic_symbol_lookup = dict((k, v.lower()) for (k, v) in _aliphatic_symbol_lookup.items())
_aliphatic_symbol_lookup["#0"] = "*"

_smarts_atom_pat = re.compile(r"""
\[
([^;]+)  # 1 = atom symbol
;
X(\d+)   # 2 = X
;
H(\d+)   # 3 = number of hydrogens
;
([+-]\d+)   # 4 = charge
;
!?R
(:\d+)?    # 5 = attachment point
\]
(([0-9]|%[0-9][0-9])*)  # 6 = any possible closures
""", re.X)
def get_environment_pseudosmiles_from_smarts(smarts):
    matches = list(_smarts_atom_pat.finditer(smarts))
    matches.reverse()
    
    pat = Chem.MolFromSmarts(smarts)
    assert pat.GetNumAtoms() == len(matches), (smarts, matches)
    
    # Work in reverse order to be able to modify the SMILES string in-place-ish.
    smiles = smarts
    for match, pat_atom in zip(matches, reversed(pat.GetAtoms())):
        bonds = pat_atom.GetBonds()

        hcount = int(match.group(3))
        
        extra_count = int(match.group(2)) - hcount - len(bonds)
        
        element_term = match.group(1)
        if element_term[0] == "#":
            # Tried but failed earlier to figure out aromaticity.
            # Try again.
            # If any aromatic bonds were included in the SMARTS
            # then the atom must be aromatic
            if any(bond.GetSmarts() == ":" for bond in bonds):
                element_term = _aromatic_symbol_lookup[element_term]
            elif extra_count < 2:
                # This must be aliphatic.
                # Aromatic atoms must have at least two aromatic bonds.
                element_term = _aliphatic_symbol_lookup[element_term]
            else:
                # Still unknown! What should I do?
                pass
            
        else:
            # Leave it as-is
            pass

        if hcount == 0:
            hydrogens = ""
        elif hcount == 1:
            hydrogens = "H"
        else:
            hydrogens = "H" + str(hcount)
            
        charge = match.group(4)
        if charge == "+0":
            charge = ""

        attachment_label = match.group(5)
        if attachment_label is None:
            attachment_label = ""
        closures = match.group(6)
        smiles_term = "[" + element_term + hydrogens + charge + attachment_label + "]" + closures + "(~*)"*extra_count
        
        smiles = smiles[:match.start()] + smiles_term + smiles[match.end():]
    
    return smiles


# The input is a SMILES for the constant term, like
#   CCC([*:2])CCC([*:1])CCCCCC[*:3]
# The centers are the attachment points 1, 2, and 3.
#


def compute_constant_environment_from_centers(centers, min_radius=0, max_radius=5):
    env_fps = []
    smarts_list = get_environment_smarts_list_for_center(centers, max_radius=max_radius)
    parent_smarts = None
    for radius in range(min_radius, max_radius+1):
        env_smarts = smarts_list[radius]
        env_smi = get_environment_pseudosmiles_from_smarts(env_smarts)
        env_fps.append(
            EnvironmentFingerprint(radius, env_smarts, env_smi, parent_smarts)
            )
        parent_smarts = env_smarts

    return env_fps

def _add_labels(constant_smiles):
    # If the constant is unlabeled, get the labeled version.
    if "*" in constant_smiles:
        # Handle differences in * representation in pre/post 2018 versions of RDKit
        if "[*]" in constant_smiles:
            wildcard = "[*]"
            width = 3
        else:
            wildcard = "*"
            width = 1
        # The conversion is direct; it's already canonical, in the order 1, 2, 3.
        # The following works for 1, 2, or 3 attachment points. If n=1 then
        # the last two replace() functions do nothing.
        # The main problem that if the wildcard is "*" then I can't do a
        # simple string substitution otherwise I'll match the previous *.
        i = 0
        i = constant_smiles.find(wildcard, i)
        if i >= 0:
            constant_smiles = constant_smiles[:i] + "[*:1]" + constant_smiles[i + width :]
            i = constant_smiles.find(wildcard, i + 5)
            if i >= 0:
                constant_smiles = constant_smiles[:i] + "[*:2]" + constant_smiles[i + width :]
                i = constant_smiles.find(wildcard, i + 5)
                if i >= 0:
                    constant_smiles = constant_smiles[:i] + "[*:3]" + constant_smiles[i + width :]

    return constant_smiles

def compute_constant_center_smarts_list(constant_smiles, min_radius=0, max_radius=5):
    orig_constant_smiles = constant_smiles
    constant_smiles = _add_labels(constant_smiles)
    env_centers = find_centers(constant_smiles)

    all_center_smarts_list = []
    seen = set()
    for smarts in get_environment_smarts_list_for_center(env_centers, min_radius=min_radius, max_radius=max_radius):
        if smarts in seen:
            continue
        seen.add(smarts)
        all_center_smarts_list.append(smarts)

    return all_center_smarts_list


######

# The "reorder" describes the changes in the the core. For example:

_invert_order_table = {
    None: None,
    "1": None,
    "12": None,
    "21": (1, 0),
    "123": None,
    "132": (0, 2, 1),
    "213": (1, 0, 2),
    "231": (2, 0, 1),
    "312": (1, 2, 0),
    "321": (2, 1, 0),
}



def _make_smarts(*center_smarts_list):
    N = len(center_smarts_list)
    if N == 1:
        return center_smarts_list[0]

    if N == 2:
        A, B = center_smarts_list
        A = A.replace(":2", ":1")
        B = B.replace(":1", ":2")
        return A + "." + B
    
    if N == 3:
        A, B, C = center_smarts_list
        A = A.replace(":2", ":1").replace(":3", ":1")
        B = B.replace(":1", ":2").replace(":3", ":2")
        C = C.replace(":1", ":3").replace(":2", ":3")
        return A + "." + B + "." + C
    
    raise AssertionError(("too many terms", center_smarts_list))

def compute_possible_smarts_environments(center_smarts_list, symmetry_class, reorder=None):
    # This is complex.
    # The center_fps are based on the constant/context fragments, which have the
    # attachment orders of 1, 2, 3.

    # The attachment orders of the core SMILES is arbitrary. Not only
    # might the (semi-)canonical SMILES be  [*:2]C([*:3])N[*:1],
    # but the attachments might be reorderd into one of the n! forms,
    # like [*:1]C([*:3])N[*:2] or [*:2]C([*:1])N[*:3]. I need to
    # unorder the permtuations to work with the 'true' numbering.

    # One I have that, I need to find all the ways to match the
    # core attachents the context. Consider:
    #   [*:2]CN([*:3])N[*:1]
    # This has a symmetry class of "122" because the *:1 is in its
    # own symmetry group and the *:2 and *:3 are in the same group.
    # (The symmetry group labels are in attachment number order, not
    # position in the SMILES.)

    # If this is attached to "[:*1]N.[:*2]N.[:*3]P" then there's
    # some ambiguity. The environment context should be the same
    # as if it were attached to "[:*1]N.[:*2]P.[:*3]N", because the
    # :*2 and *:3 are symmetric.

    # This means I can have to enumerate all of the possiblities, to
    # find all of the possible environments.

    invert_order = _invert_order_table[reorder]
    if invert_order is not None:
        center_smarts_list = [center_smarts_list[i] for i in invert_order]

    smarts_list = set()
    if symmetry_class in ("1", "12", "123"):
        smarts_list.add(_make_smarts(*center_smarts_list))
    elif symmetry_class == "11":
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[1]))
        smarts_list.add(_make_smarts(center_smarts_list[1], center_smarts_list[0]))

    elif symmetry_class == "111":
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[1], center_smarts_list[2]))
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[2], center_smarts_list[1]))
        smarts_list.add(_make_smarts(center_smarts_list[1], center_smarts_list[0], center_smarts_list[2]))
        smarts_list.add(_make_smarts(center_smarts_list[1], center_smarts_list[2], center_smarts_list[0]))
        smarts_list.add(_make_smarts(center_smarts_list[2], center_smarts_list[0], center_smarts_list[1]))
        smarts_list.add(_make_smarts(center_smarts_list[2], center_smarts_list[1], center_smarts_list[0]))

    elif symmetry_class == "112":
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[1], center_smarts_list[2]))
        smarts_list.add(_make_smarts(center_smarts_list[1], center_smarts_list[0], center_smarts_list[2]))

    elif symmetry_class == "122":
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[1], center_smarts_list[2]))
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[2], center_smarts_list[1]))

    elif symmetry_class == "121":
        smarts_list.add(_make_smarts(center_smarts_list[0], center_smarts_list[1], center_smarts_list[2]))
        smarts_list.add(_make_smarts(center_smarts_list[2], center_smarts_list[1], center_smarts_list[0]))

    else:
        raise AssertionError("I forgot one: %r" % (symmetry_class,))

    return list(smarts_list)

def get_all_possible_smarts(all_center_smarts_list, symmetry_class, permutation):
    possible_envs = set()
    for center_smarts in all_center_smarts_list:
        possible_envs.update(compute_possible_smarts_environments(
            center_smarts.split("."), symmetry_class, permutation))
    return sorted(possible_envs)
