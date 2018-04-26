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

# Determine the environment around the attachment point(s) of the constant
# part using circular fingerprints.

import sys
import hashlib
import binascii

from rdkit import Chem
from rdkit.Chem import AllChem

from collections import namedtuple

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
    __slots__ = ("radius", "fingerprint")
    def __init__(self, radius, fingerprint):
        self.radius = radius
        self.fingerprint = fingerprint
    

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

        # Don't need to clear the property because RDKit's Morgan fingerprint
        # doesn't use it.
        #atom.ClearProp("molAtomMapNumber")
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
            ordered_centers = centers[1], # this is a one element tuple 
        except KeyError:
            raise ValueError("Expecting attachment point *:1 in context SMILES %r" % (smiles,))

    return EnvironmentCenters(mol, ordered_centers)

# Get each of the atoms counts radii as it's being computed.
def iter_num_atoms_for_radii(centers, max_radius):
    return _iter_num_atoms_for_radii(centers.mol, max_radius, centers.atom_ids)

# My thought is to use this for testing.
def _iter_num_atoms_for_radii(mol, max_radius, start_atoms): 
    atom_objs = list(mol.GetAtoms())
    
    unique_atoms = set(start_atoms)
    assert len(start_atoms) == len(unique_atoms), "duplicate start atom"
    ignore_atoms = set(a for a in start_atoms if not is_heavy_atom(mol.GetAtomWithIdx(a)))
    
    yield (len(unique_atoms) - len(ignore_atoms))
    
    border_atoms = unique_atoms.copy()
    
    for radius in range(max_radius):  # up to and including max_radius
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

def find_center_fingerprints(centers, radius):
    center_fps = []
    for atom_idx in centers.atom_ids:
        fp = AllChem.GetMorganFingerprint(centers.mol, radius, fromAtoms=[atom_idx])
        # Hash so everything so they are a constant 32 binary bytes.
        center_fp = hashlib.sha256(fp.ToBinary()).digest()
        center_fps.append(center_fp)
    return center_fps

def _make_fp(*center_fps):
    concat_fp = hashlib.new("sha256")
    for center_fp in center_fps:
        concat_fp.update(center_fp)
        
    sha2_digest = concat_fp.digest()
    env_fp = binascii.b2a_base64(sha2_digest)
    assert env_fp[-2:] == b"=\n", env_fp
    return fingerprint_to_text(env_fp[:-2]) # trim the "=\n"


def find_environment_fingerprint(centers, radius):
    # I want a unique fingerpint for the given environment.

    # The attachment points are in canonical order so I can
    # concatenate the circular fingerprints, making sure that there's
    # only one way to produce a given concatenation. (There's a low
    # chance that simple concatenation of two ToBinary() strings,
    # which are variable length, might be ambiguous with a
    # concatenation of other strings.)

    # Fix the ambiguity problem by using the SHA2 of the ToBinary() so
    # all circular fingerprints are a constant length.

    # The result is far longer than it needs to be. I only need to
    # test two environments for identify. I don't need to know the
    # per-center breakdown. Use sha2 again to solve that problem.

    # TODO: the fingerprints at r=0 are constant based on the number of cuts.
    # Would it be faster to have special code for that case?
    concat_fp = _make_fp(*find_center_fingerprints(centers, radius))
    return concat_fp

# The input is a SMILES for the constant term, like
#   CCC([*:2])CCC([*:1])CCCCCC[*:3]
# The centers are the attachment points 1, 2, and 3.
# 


def compute_constant_environment_from_centers(centers, max_radius=5):
    env_fps = []
    for radius in range(max_radius+1):
        fingerprint = find_environment_fingerprint(centers, radius)
        env_fps.append(
            EnvironmentFingerprint(radius, fingerprint)
            )

    return env_fps

def compute_constant_center_fingerprints(constant_smiles, min_radius=0, max_radius=5):
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
            constant_smiles = constant_smiles[:i] + "[*:1]" + constant_smiles[i+width:]
            i = constant_smiles.find(wildcard, i+5)
            if i >= 0:
                constant_smiles = constant_smiles[:i] + "[*:2]" + constant_smiles[i+width:]
                i = constant_smiles.find(wildcard, i+5)
                if i >= 0:
                    constant_smiles = constant_smiles[:i] + "[*:3]" + constant_smiles[i+width:]

    env_centers = find_centers(constant_smiles)

    all_center_fps = []
    for radius in range(min_radius, max_radius+1):
        all_center_fps.append(find_center_fingerprints(env_centers, radius))

    return all_center_fps


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


def compute_possible_environments(center_fps, symmetry_class, reorder=None):
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
        center_fps = [center_fps[i] for i in invert_order]

    fps = set()
    if symmetry_class in ("1", "12", "123"):
        fps.add(_make_fp(*center_fps))
    elif symmetry_class == "11":
        fps.add(_make_fp(center_fps[0], center_fps[1]))
        fps.add(_make_fp(center_fps[1], center_fps[0]))

    elif symmetry_class == "111":
        fps.add(_make_fp(center_fps[0], center_fps[1], center_fps[2]))
        fps.add(_make_fp(center_fps[0], center_fps[2], center_fps[1]))
        fps.add(_make_fp(center_fps[1], center_fps[0], center_fps[2]))
        fps.add(_make_fp(center_fps[1], center_fps[2], center_fps[0]))
        fps.add(_make_fp(center_fps[2], center_fps[0], center_fps[1]))
        fps.add(_make_fp(center_fps[2], center_fps[1], center_fps[0]))

    elif symmetry_class == "112":
        fps.add(_make_fp(center_fps[0], center_fps[1], center_fps[2]))
        fps.add(_make_fp(center_fps[1], center_fps[0], center_fps[2]))

    elif symmetry_class == "122":
        fps.add(_make_fp(center_fps[0], center_fps[1], center_fps[2]))
        fps.add(_make_fp(center_fps[0], center_fps[2], center_fps[1]))

    elif symmetry_class == "121":
        fps.add(_make_fp(center_fps[0], center_fps[1], center_fps[2]))
        fps.add(_make_fp(center_fps[2], center_fps[1], center_fps[0]))

    else:
        raise AssertionError("I forgot one: %r" % (symmetry_class,))

    return list(fps)

## def generate_all_possible_constant_fingerprints(constant_smiles):
##     # The n circular fingerprints around each of the n attachment points, for each possible radius
##     all_center_fps = compute_constant_center_fingerprints(constant_smiles)
##     n = constant_smiles.count("*")
##     symclass = "1" * n
##     fingerprints = set()
##     for center_fps in all_center_fps:
##         fingerprints.update(compute_possible_environments(center_fps, symclass))
##     return fingerprints

def get_all_possible_fingerprints(all_center_fps, symmetry_class, permutation):
    possible_envs = set()
    for center_fps in all_center_fps:
        possible_envs.update(compute_possible_environments(
            center_fps, symmetry_class, permutation))
    return sorted(possible_envs)
