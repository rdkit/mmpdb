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


class CutSmarts(object):
    def __init__(self, name, smarts, description):
        self.name = name
        self.smarts = smarts
        self.description = description


cut_smarts_aliases_by_name = {}

cut_smarts_aliases = [
        CutSmarts(
            "default",
            "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]",
            "Cut all C-[!H] non-ring single bonds except for Amides/Esters/Amidines/Sulfonamides "
            "and CH2-CH2 and CH2-CH3 bonds"),
     
        CutSmarts(
            "cut_AlkylChains",
            "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1]",
            "As default, but also cuts CH2-CH2 and CH2-CH3 bonds"),
            
        CutSmarts(
            "cut_Amides",
            "[#6+0]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]",
            "As default, but also cuts [O,N]=C-[O,N] single bonds"),

        CutSmarts(
            "cut_all",
            "[#6+0]!@!=!#[!#0;!#1]",
            "Cuts all Carbon-[!H] single non-ring bonds. Use carefully, this will create a lot of cuts"),

        CutSmarts(
            "exocyclic",
            "[R]!@!=!#[!#0;!#1]",
            "Cuts all exocyclic single bonds"),

        CutSmarts(
            "exocyclic_NoMethyl",
            "[R]!@!=!#[!#0;!#1;!$([CH3])]",
            "Cuts all exocyclic single bonds apart from those connecting to CH3 groups"),
    ]


for alias in cut_smarts_aliases:
    cut_smarts_aliases_by_name[alias.name] = alias


def get_epilog(option_name, aliases):
    lines = ["The " + option_name + " argument supports the following short-hand aliases:"]
    for alias in aliases:
        lines.append("  '%s': %s" % (alias.name, alias.description))
        lines.append("     smarts: %s" % (alias.smarts,))
    return "\n".join(lines)+"\n"
