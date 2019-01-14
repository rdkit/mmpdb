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


class FragmentRecord(object):
    __slots__ = ("id", "input_smiles", "num_normalized_heavies", "normalized_smiles", "fragments")
    errmsg = None

    def __init__(self, id, input_smiles, num_normalized_heavies, normalized_smiles, fragments):
        self.id = id
        self.input_smiles = input_smiles
        self.num_normalized_heavies = num_normalized_heavies
        self.normalized_smiles = normalized_smiles
        self.fragments = fragments

    def __repr__(self):
        return "FragmentRecord(%r, %r, %d, %r, %r)" % (
            self.id, self.input_smiles, self.num_normalized_heavies, self.normalized_smiles,
            self.fragments)


class FragmentErrorRecord(object):
    __slots__ = ("id", "input_smiles", "errmsg", "fragments")
    num_normalized_heavies = normalized_mol = None
    
    def __init__(self, id, input_smiles, errmsg):
        self.id = id
        self.input_smiles = input_smiles
        self.errmsg = errmsg
        self.fragments = []  # Don't share a mutable list with other instances

    def __repr__(self):
        return "FragmentErrorRecord(%r, %r, %r)" % (self.id, self.input_smiles, self.errmsg)


class Fragment(object):
    __slots__ = ("num_cuts", "variable_symmetry_class", "num_variable_heavies", "variable_smiles",
                 "num_constant_heavies", "constant_smiles", "constant_with_H_smiles")

    def __init__(self, num_cuts,
                 variable_symmetry_class, num_variable_heavies, variable_smiles,
                 num_constant_heavies, constant_smiles, constant_with_H_smiles):
        ## assert num_cuts in (1, 2, 3)
        ## assert len(variable_symmetry_class) in (1,2,3), variable_symmetry_class
        ## assert num_variable_heavies > 0
        ## assert variable_smiles
        ## assert num_constant_heavies > 0
        ## assert constant_smiles

        ## if "[*:2]" in constant_smiles:
        ##     assert constant_smiles.index("[*:1]") < constant_smiles.index("[*:2]"), constant_smiles
        ##     if "[*:3]" in constant_smiles:
        ##         assert constant_smiles.index("[*:2]") < constant_smiles.index("[*:3]"), constant_smiles
        
        self.num_cuts = num_cuts
        self.variable_symmetry_class = variable_symmetry_class
        self.num_variable_heavies = num_variable_heavies
        self.variable_smiles = variable_smiles
        self.num_constant_heavies = num_constant_heavies
        self.constant_smiles = constant_smiles
        self.constant_with_H_smiles = constant_with_H_smiles
        
    def __repr__(self):
        return "Fragment(%d, %r, %d, %r, %d, %r, %r)" % (
            self.num_cuts, self.variable_symmetry_class,
            self.num_variable_heavies, self.variable_smiles,
            self.num_constant_heavies, self.constant_smiles, self.constant_with_H_smiles)
    
        
class FragmentValueError(ValueError):
    def __init__(self, name, value, reason):
        self.name = name
        self.value = value
        self.reason = reason

    def __str__(self):
        return "Error with %s (%r): %s" % (self.name, self.value, self.reason)

    def __repr__(self):
        return "FragmentValueError(%r, %r, %r)" % (
            self.name, self.value, self.reason)


class FragmentFormatError(ValueError):
    def __init__(self, reason, location):
        self.reason = reason
        self.location = location

    def __repr__(self):
        return "FragmentFormatError(%r, %r)" % (
            self.reason, self.location)
    
    def __str__(self):
        return "%s, %s" % (self.reason, self.location.where())
