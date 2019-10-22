# Convert a molecule fragment into a SMARTS pattern.
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
import sys
from rdkit import Chem
import argparse

from . import command_support

# Match the atom terms. These are all in brackets.
_atom_term = re.compile(r"\[([^]]+)\]")

def fragment_mol_to_smarts(mol):
    """mol -> SMARTS matching the fragment R-group

    'mol' contain a single wildcard atom ("*") with one
    single bond to the rest of the R-group.

    The molecule will be modified. Pass in a copy if
    you do not want the modification.
    """

    wildcard_idx = -1
    suffixes = []

    for idx, atom in enumerate(mol.GetAtoms()):
        if atom.HasProp("molAtomMapNumber"):
            atom_map = atom.GetProp("molAtomMapNumber")
            if atom_map != "1":
                raise ValueError("atom maps are not supported (atom %d has atommap %r)" % (
                    idx, atom_map))
        
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
    converted_smi = Chem.MolToSmiles(mol,
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
        super(ParseError, self).__init__(msg)
        self.location = location
        self._where = location.where()
    def __str__(self):
        return "%s at %s" % (self.msg, self._where)

class ConversionError(ValueError):
    def __init__(self, msg, location, extra=None):
        super(ConversionError, self).__init__(msg)
        self.location = location
        self.extra = None
        self._where = location.where()
    def __str__(self):
        if self.extra is None:
            return "%s at %s" % (self.msg, self._where)
        else:
            return "%s at %s: %s" % (self.msg, self._where, self.extra)

class Record(object):
    def __init__(self, smiles, id=None):
        self.smiles = smiles
        self.id = id
    def __repr__(self):
        return "Record(%r, id=%r)" % (self.smiles, self.id)
        
def parse_fragment_file(infile, location=None):
    if location is None:
        location = FileLocation(getattr(infile, "name", "<unknown>"))

    recno = 0
    for lineno, line in enumerate(infile, 1):
        location.lineno = lineno
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
        recno += 1
        location.recno = recno
        yield rec
            
######

class Location(object):
    def __init__(self, source, record_id=None):
        self.source = source
        self.record_id = record_id
        self.recno = 0
        
    def where(self):
        if self.record_id is None:
            return "%s" % (self.source,)
        else:
            return "%s, record %r" % (self.source, self.record_id)

class FileLocation(Location):
    def __init__(self, source):
        super(Location, self).__init__(source)
        self.lineno = 0

    def where(self):
        if self.record_id is None:
            return "%s, line %d" % (self.source, self.lineno)
        else:
            return "%s, line %d, record %r" % (self.source, self.lineno, self.record_id)

    def __repr__(self):
        return "FileLocation(%r, record_id=%r, lineno=%d)" % (
            self.source, self.record_id, self.lineno)

class ListLocation(Location):
    def __init__(self, source, initial_index=0):
        super(ListLocation, self).__init__(source)
        self.index = initial_index
        
    def where(self):
        if self.record_id is None:
            return "%s #%s" % (self.source, self.index)
        else:
            return "%s #%s, record %r" % (self.source, self.index, self.record_id)

    def __repr__(self):
        return "ListLocation(%r, record_id=%r, index=%d)" % (
            self.source, self.record_id, self.index)

def iter_smiles_list(smiles_iter, location):
    for i, smiles in enumerate(smiles_iter, location.index):
        location.index = location.recno = i
        yield Record(smiles, None)


        
def iter_smiles_as_smarts(record_reader, location, explain=None, all_mols=None):
    if explain is None:
        explain = command_support.no_explain
        
    for record in record_reader:
        location.record_id = record.id
        smiles = record.smiles
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ConversionError(
                "Cannot parse SMILES (%r)" % (smiles,),
                location
                )
        try:
            smarts = fragment_mol_to_smarts(mol)
        except ValueError as err:
            raise ConversionError(
                "Cannot convert SMILES (%r)" % (smiles,),
                location_error,
                str(err))
        
        explain("#%d: converted SMILES %r to SMARTS %r" % (location.recno, smiles, smarts))
        if all_mols is not None:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                raise ConversionError("SMARTS failure for %r" % (smarts,),
                                          location,
                                          "using SMILES %r" % (smiles,))
            
            if not mol.HasSubstructMatch(pat):
                raise ConversionError("SMARTS %r does not match SMILES %r" % (smarts, smiles),
                                          location)
            all_mols.append((mol, location.where(), smiles))
            explain("#%d passed the self-check" % (location.recno,))
            
        yield smarts

        
###### Command-line code

def die(msg):
    sys.stdout.write(msg + "\n")
    raise SystemExit(1)

def frag2smarts_command(parser, args):
    check = args.check
    explain = command_support.get_explain(args.explain)
    
    filename = "<unknown>"
    close = None
    
    if args.smiles is not None:
        if args.fragment_filename is not None:
            parser.error("Cannot specify both a fragment filename and a fragment --smiles")
        location = ListLocation("command-line SMILES", 1)
        explain("Using SMILES from the command-line")
        record_reader = iter_smiles_list(args.smiles, location)
        
    elif args.fragment_filename is not None:
        filename = args.fragment_filename
        explain("Reading SMILES from %r" % (filename,))
        location = FileLocation(args.input)
        try:
            f = open(args.input)
        except OSError as err:
            die("Cannot open input file: %s" % (err,))
        close = f.close
        record_reader = parse_fragment_file(f, location)
        
    else:
        explain("Reading SMILES from <stdin>")
        location = FileLocation("<stdin>")
        record_reader = parse_fragment_file(sys.stdin, location)

    if check:
        all_mols = []
    else:
        all_mols = None
        
    outfile = sys.stdout

    iter_smarts = iter_smiles_as_smarts(record_reader, location, explain, all_mols)

    all_smarts = None
    
    try:
        if args.single:
            for smarts in iter_smarts:
                outfile.write(smarts + "\n")
        else:
            all_smarts = []
            for smarts in iter_smarts:
                assert smarts.startswith("*-!@"), (smarts, location)
                all_smarts.append(smarts[4:])
            if not all_smarts:
                die("Cannot make a SMARTS: no SMILES strings found in %s" % (location.source,))
                
    except ParseError as err:
        die("Cannot parse input file: %s" % (err,))
    except ConversionError as err:
        die(str(err))
    finally:
        if close is not None:
            close()


    if not args.single:
        # Create the recursive SMILES
        rest = "[" + ",".join("$(%s)" % smarts for smarts in all_smarts) + "]"
        smarts = "*-!@" + rest
        
        try:
            if check:
                explain("Checking that the SMARTS matches all of the input molecules")
                all_pat = Chem.MolFromSmarts(smarts)
                if all_pat is None:
                    die("Cannot process final SMARTS: %r" % (smarts,))
                
                for i, (mol, where, smiles) in enumerate(all_mols):
                    if not mol.HasSubstructMatch(all_pat):
                        die("final SMARTS does not match SMILES from %s (%r)" % (
                            where, smiles))
                    explain("checked #%d" % (i,))
        finally:
            outfile.write(smarts + "\n")
            
    outfile.flush()
    

if __name__ == "__main__":
    main()
