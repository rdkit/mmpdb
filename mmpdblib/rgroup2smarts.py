# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2019, Andrew Dalke Scientific, AB
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
#    * Neither the name of Andrew Dalke Scientific. nor the names of
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

# Convert a R-group SMILES into a SMARTS pattern.
# Must be rooted at "*", with one single bond.

# The algorithm is:
#  - figure out the valence and hydrogen count for each atom
#  - convert the molecule into an isomeric SMILES, with
#      hydrogens on the atom terms
#  - convert the SMILES into a SMARTS:
#     - use a regexp to find the atom terms
#     - insert with the valence

# For examples:
#   *-Cl -> *!@[Clv1]
#   *-O  -> *!@[OH1v2]

import re
from rdkit import Chem

from . import fileio

# Match the atom terms. These are all in brackets.
_atom_term = re.compile(r"\[([^]]+)\]")


def rgroup_mol_to_smarts(mol):
    """mol -> SMARTS matching the R-group

    'mol' contain a single wildcard atom ("*") with one
    single bond to the rest of the R-group.

    The molecule will be modified. Pass in a copy if
    you do not want the modification.
    """

    wildcard_idx = -1
    suffixes = []

    n = len(Chem.GetMolFrags(mol))
    if n > 1:
        raise ValueError("more than one fragment found")
    # Could also check for n == 0. Decided to leave that
    # for the check for no wildcard atom.

    for idx, atom in enumerate(mol.GetAtoms()):
        if atom.HasProp("molAtomMapNumber"):
            atom_map = atom.GetProp("molAtomMapNumber")
            raise ValueError(f"atom maps are not supported (atom {idx} has atom map {atom_map!r})")

        # Figure out which atom is the wildcard
        atomno = atom.GetAtomicNum()
        if atomno == 0:
            if wildcard_idx == -1:
                wildcard_idx = idx
            else:
                raise ValueError("more than one wildcard atom")
            bonds = list(atom.GetBonds())
            if not bonds:
                raise ValueError("wildcard atom not bonded to anything")
            if len(bonds) != 1:
                raise ValueError("wildcard atom must only have one bond")
            if bonds[0].GetBondType() != 1:
                raise ValueError("wildcard atom not bonded via a single bond")

            if atom.GetTotalNumHs():
                raise ValueError("wildcard atom must not have implicit hydrogens")
            if atom.GetFormalCharge():
                raise ValueError("wildcard atom must be uncharged")
            # I could check for the chiral flag, but RDKit won't use it.

        # The suffix is after the atom term.
        v = atom.GetTotalValence()
        suffix = "v" + str(v)

        # If an atom has no hydrogens then MolToSmiles(allHsExplicit=True)
        # will not include an H term. When used as a SMARTS, this means
        # the atom term will match anything, when I want to to match nothing.
        # I need to also include an H0 in the SMARTS.
        if not atom.GetTotalNumHs():
            suffix = "H0" + suffix
        suffixes.append(suffix)

    if wildcard_idx == -1:
        raise ValueError("no wildcard atom found")

    # allHsExplicit=True so all of the atom terms are in [brackets].
    # allBondsExplicit=True to distinguish between aromatic and single bonds
    #    (which SMILES handles as an implicit bond)
    converted_smi = Chem.MolToSmiles(
        mol,
        allBondsExplicit=True,
        allHsExplicit=True,
        rootedAtAtom=wildcard_idx,
    )

    # Figure out how to get from the canonical Kekule SMILES output order
    # back to the moleule order.
    # The MolToSmiles() set the '_smilesAtomOutputOrder'property, which
    # is a string like '[0,1,2,3,4,5,6,7,8,13,9,10,11,12,]' which maps
    # between the atom order in the molecule and the output order
    output_order_str = mol.GetProp("_smilesAtomOutputOrder")
    assert output_order_str[0] == "[" and output_order_str[-2:] == ",]", output_order_str
    output_order_map = map(int, output_order_str[1:-2].split(","))
    invert_order = dict(enumerate(output_order_map))

    ## print(converted_smi)
    ## print(suffixes)

    # Use a regular expression to

    # This emulates 'nonlocal' support so I can get
    # the position of the re.sub() index even under
    # Python 2.7.
    sub_index_nonlocal = [0]

    def rename_atoms(m):
        sub_index = sub_index_nonlocal[0]
        sub_index_nonlocal[0] = sub_index + 1

        original_index = invert_order[sub_index]
        suffix = suffixes[original_index]
        g = m.group(1)
        return "[" + g + suffix + "]"

    output_smi = _atom_term.sub(rename_atoms, converted_smi)
    assert output_smi.startswith("[*H0v1]-"), output_smi
    return "*-!@" + output_smi[8:]


class ParseError(ValueError):
    def __init__(self, msg, location):
        super(ParseError, self).__init__(msg, location)
        self.msg = msg
        self.location = location
        self._where = location.where()

    def __str__(self):
        return f"{self.msg} at {self._where}"


class ConversionError(ValueError):
    def __init__(self, msg, location, extra=None):
        super(ConversionError, self).__init__(msg, location, extra)
        self.msg = msg
        self.location = location
        self.extra = extra
        self._where = location.where()

    def __str__(self):
        if self.extra is None:
            return f"{self.msg} at {self._where}"
        else:
            return f"{self.msg} at {self._where}: {self.extra}"


class Record(object):
    def __init__(self, smiles, id=None):
        self.smiles = smiles
        self.id = id

    def __repr__(self):
        return "Record({self.smiles!r}, id={self.id!r})"


def parse_rgroup_file(infile, location=None):
    if location is None:
        location = FileLocation(getattr(infile, "name", "<unknown>"))
    recno = 0
    lineno = 0

    def get_recno():
        return lineno

    def get_lineno():
        return lineno

    location.register(get_recno=get_recno, get_lineno=get_lineno)

    try:
        for lineno, line in enumerate(infile, 1):
            if line == "\n":
                raise ParseError("no SMILES found", location)

            if line[:1] in "\r\v\t ":
                raise ParseError("expected SMILES at start of line", location)

            terms = line.split(None, 1)
            if not terms:
                raise ParseError("no SMILES found", location)
            elif len(terms) == 1:
                rec = Record(terms[0], None)
            else:
                rec = Record(terms[0], terms[1].rstrip("\n\r"))
            yield rec
    finally:
        location.save(recno=recno, lineno=lineno)


######


class FileLocation(fileio.Location):
    def where(self):
        msg = f"{self.filename!r}, line {self.lineno}"
        if self.record_id is not None:
            msg += f", record {self.record_id!r}"
        return msg


class ListLocation(fileio.Location):
    def where(self):
        if self.record_id is None:
            return f"{self.filename} #{self.recno}"
        else:
            return f"{self.source} #{self.recno}, record {self.record_id!r}"


def iter_smiles_list(smiles_iter, location):
    recno = location.recno

    def get_recno():
        return recno

    location.register(get_recno=get_recno)
    try:
        for recno, smiles in enumerate(smiles_iter, recno):
            yield Record(smiles, None)
    finally:
        location.save(recno=recno)


def iter_smiles_as_smarts(record_reader, location, explain=None, all_mols=None):
    if explain is None:

        def explain(msg, *args):
            pass

    for record in record_reader:
        location.record_id = record.id
        smiles = record.smiles
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ConversionError(f"Cannot parse SMILES ({smiles!r})", location)
        try:
            smarts = rgroup_mol_to_smarts(mol)
        except ValueError as err:
            raise ConversionError(f"Cannot convert SMILES ({smiles!r})", location, str(err))

        explain(f"#{location.recno}: converted SMILES {smiles!r} to SMARTS {smarts!r}")
        if all_mols is not None:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                raise ConversionError(
                    f"SMARTS failure for {smarts!r}",
                    location,
                    "using SMILES {smiles!r}",
                )

            if not mol.HasSubstructMatch(pat):
                raise ConversionError(
                    "SMARTS {smarts!r} does not match SMILES {smiles!r}",
                    location,
                )
            all_mols.append((mol, location.where(), smiles))
            explain(f"#{location.recno} passed the self-check")

        yield smarts


def make_recursive_smarts(smarts_list):
    terms = []
    for smarts in smarts_list:
        if not smarts.startswith("*-!@"):
            raise ValueError("invalid prefix: {smarts!r}")

        terms.append("$(" + smarts[4:] + ")")
    return "*-!@[" + ",".join(terms) + "]"


def get_recursive_smarts_from_cut_rgroups(rgroups, source="rgroup", offset=0):
    location = ListLocation(source)
    location.save(recno=offset)
    record_reader = iter_smiles_list(rgroups, location)
    iter_smarts = iter_smiles_as_smarts(record_reader, location)
    return make_recursive_smarts(iter_smarts)


def get_recursive_smarts_from_cut_filename(filename):
    location = FileLocation(filename)
    with open(filename) as infile:
        record_reader = parse_rgroup_file(infile, location)
        iter_smarts = iter_smiles_as_smarts(record_reader, location)
        return make_recursive_smarts(iter_smarts)
