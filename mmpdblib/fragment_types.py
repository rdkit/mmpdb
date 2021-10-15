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
# 3) fragment_io.py -- defines the mapping from SQL to the data types

# I tried using SQLAlchemy to merge them into one but that
# added a ~30% overhead to fragment generation.

import dataclasses
from dataclasses import dataclass
from typing import List, Optional, Literal

@dataclass
class FragmentOptions:
    max_heavies: Optional[int]
    max_rotatable_bonds: Optional[int]
    rotatable_smarts: str
    cut_smarts: str

    num_cuts: Literal[1,2,3]
    method: Literal["chiral"]
    salt_remover: str
    min_heavies_per_const_frag: int

    def to_dict(self):
        return dataclasses.asdict(self)
    
    # For backwards compatability (still needed?) -- TODO: remove
    def to_text_settings(self):
        def _none(x):
            return "none" if x is None else str(x)
        return tuple((name, _none(value))
                         for (name, value) in dataclasses.asdict(self).items())

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
        return "%s.%s.%s" % (self.attachment_order, self.variable_smiles, self.constant_smiles)
    
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
