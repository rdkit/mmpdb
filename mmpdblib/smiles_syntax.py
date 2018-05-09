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

# I call this "smiles_syntax.py"
from __future__ import print_function

import re

# Match a '*' in the different forms that might occur,
# including with directional single bonds inside of ()s.
_wildcard_regex = " |\n".join(re.escape(regex) for regex in
  ("*", "[*]", "(*)", "([*])", "(/*)", "(/[*])", "/*", "/[*]", "(\\*)", "(\\[*])", "\\*", "\\[*]"))
_wildcard_pattern = re.compile(_wildcard_regex, re.X)

# Match the SMILES for an atom, followed by its closures
_atom_pattern = re.compile(r"""
(
 Cl? |             # Cl and Br are part of the organic subset
 Br? |
 [NOSPFIbcnosp*] |  # as are these single-letter elements
 \[[^]]*\]         # everything else must be in []s
)
""", re.X)

def convert_wildcards_to_closures(smiles, offsets=None):
    # This is designed for RDKit's canonical SMILES output. It does
    # not handle all possible SMILES inputs.
    if offsets is None:
        # Use 0, 1, 2, ... up to the number of '*'s
        offsets = range(smiles.count("*"))
    closure_terms = []
    for offset in offsets:
        if not (0 <= offset <= 9):
            raise ValueError("offset %d out of range (must be from 0 to 9)"
                             % (offset,))
        closure_terms.append("%%%02d" % (90 + offset))
    
    new_smiles = smiles
    while 1:
        # Find the first '*'. If none are left, stop.
        wildcard_match = _wildcard_pattern.search(new_smiles)
        if wildcard_match is None:
            break

        closure_term = closure_terms.pop(0)
        
        wildcard_start = wildcard_match.start()
        if wildcard_start == 0 or new_smiles[wildcard_start-1] == ".":
            # At the start of the molecule or after a ".". Need to
            # put the closure after the second atom. Find the second
            # atom. Since we only ever break on single non-ring bonds,
            # and since the first atom is a terminal atom, the second
            # atom must either be immediately after the first atom, or
            # there is a directional bond between them.
            wildcard_end = wildcard_match.end()
            second_atom_match = _atom_pattern.match(new_smiles, wildcard_end)
            if second_atom_match is None:
                # There was no atom. Is it a "/" or "\"? If so,
                # we'll need to swap the direction when we move
                # to a closure after the second atom.
                bond_dir = new_smiles[wildcard_end:wildcard_end+1]
                if bond_dir == "/":
                    direction = "\\"
                elif bond_dir == "\\":
                    direction = "/"
                else:
                    raise AssertionError(new_smiles)
                # Look for the second atom, which must exist
                second_atom_match = _atom_pattern.match(new_smiles, wildcard_end+1)
                if second_atom_match is None:
                    raise AssertionError((new_smiles, new_smiles[wildcard_end:]))
            else:
                direction = ""

            second_atom_term = second_atom_match.group(1)
            # I changed the bond configuration, so I may need to
            # invert chirality of implicit chiral hydrogens.
            if "@@H" in second_atom_term:
                second_atom_term = second_atom_term.replace("@@H", "@H")
            elif "@H" in second_atom_term:
                second_atom_term = second_atom_term.replace("@H", "@@H")

            # Reassemble the string with the wildcard term deleted and
            # the new closure inserted directly after the second atom
            # (and before any of its closures).
            new_smiles = (new_smiles[:wildcard_start] + second_atom_term
                          + direction + closure_term
                          + new_smiles[second_atom_match.end():])
                    
        else:
            # The match is somewhere inside of a molecule, so we attach
            # assign the closure to the atom it's bonded to on the left
            c = new_smiles[wildcard_start-1]
            if c == "(" or c == ")":
                # In principle, this could be something like "CCC(F)(Cl)[*]", 
                # where I would need to count the number of groups back to
                # the main atom, and flip chirality accordingly. Thankfully,
                # RDKit always puts the "[*]" terms immediately after the
                # preceeding atom, so I don't need to worry.
                raise NotImplementedError("intermediate groups not supported",
                                          new_smiles, new_smiles[wildcard_start-1:])
            
            elif c in "CNcnOS]Pos0123456789ABDEFGHIJKLMQRTUVWXYZabdefghijklmpqrtuvwxyz":
                # Double-check the the previous character looks like part of an atom.
                wildcard_term = wildcard_match.group()
                # Preserve the direction, if present
                if "/" in wildcard_term:
                    direction = "/"
                elif "\\" in wildcard_term:
                    direction = "\\"
                else:
                    direction = ""
                new_smiles = (new_smiles[:wildcard_start] + direction + closure_term
                              + new_smiles[wildcard_match.end():])

            else:
                raise AssertionError((new_smiles, c, new_smiles[wildcard_start-1:]))

    return new_smiles

##### Same thing, for labeled wildcards

_labeled_wildcard_pattern = re.compile(r"\*:([123])")


def convert_labeled_wildcards_to_closures(smiles):
    offsets = []
    def sub_function(m):
        offsets.append(int(m.group(1)))
        return "*"
    new_smiles = _labeled_wildcard_pattern.sub(sub_function, smiles)
    #print("convert_labeled_wildcards_to_closures:", smiles, new_smiles, offsets)
    return convert_wildcards_to_closures(new_smiles, offsets)


if __name__ == "__main__":
    for smiles in ("*C", "*/CO.*CN", "C*.C(*)N"):
        print(smiles, convert_wildcards_to_closures(smiles, (0,)*smiles.count("*")))
