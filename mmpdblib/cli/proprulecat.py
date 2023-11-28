"Implement the 'proprulecat' command"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2019, Relay Therapeutics, Inc.
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
#    * Neither the name of Relay Therapeutics, Inc. nor the names of
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


import click
from .click_utils import (
    command,
    die,
    GzipFile,
    add_single_database_parameters,
    positive_int,
    open_dataset_from_options_or_exit,
    get_property_names_or_error,
    )

####


import re



# Check that the given SMILES has the right number of '*'
# atoms and that they are formatted correctly.
# Returns (num_wildcard_atoms:int, is_labeled:bool)

def check_wildcards(where, smiles):
    from ..smiles_syntax import _atom_pattern
    
    num_wildcards = smiles.count("*")
    if not (1 <= num_wildcards <= 3):
        raise click.UsageError(
            f"The --{where} SMILES {smiles!r} must contain "
            f"1, 2, or 3 wildcard ('*') atoms, not {num_wildcards}"
            )

    # The only valid forms are bare '*'s, '[*]', and '[*:1]'/'[*:2]'/'[*:3]'
    atom_tokens = _atom_pattern.findall(smiles)

    # Check if they are all unlabeled, like "*" or "[*]",
    # or labeled, like "[*:1]", "[*:2]", or "[*:3]".
    
    num_unlabeled_wildcards = 0
    first_unlabeled_token = None
    seen_labeled_wildcards = set()
    for token in atom_tokens:
        # Is it an unlabeled token?
        if token == "*" or token == "[*]":
            num_unlabeled_wildcards += 1
            if first_unlabeled_token is None:
                first_unlabeled_token = token
            continue
                
        # Is it a labeled token?
        if token in ("[*:1]", "[*:2]", "[*:3]"):
            if token in seen_labeled_wildcards:
                raise click.UsageError(
                    f"Duplicate {token!r} atom found in --{where} SMILES {smiles!r}")
            
            seen_labeled_wildcards.add(token)
            continue

        # Some other token with a '*'? How odd!
        if "*" in token:
            raise click.UsageError(
                f"Wildcard atom {token!r} with unexpected format found in "
                f"--{where} SMILES {smiles!r}")
            
    # Don't mix the two types.
    if num_unlabeled_wildcards and seen_labeled_wildcards:
        first_labeled_token = sorted(seen_labeled_wildcards)[0]
        raise click.UsageError(
            f"Cannot mix both labeled ({first_labeled_token!r}) and "
            f"unlabled ({first_unlabeled_token!r}) wildcard atoms in "
            f"--{where} SMILES {smiles!r}")

    # If there are unlabled wildcards, double-check the numbers match then I'm done.
    if num_unlabeled_wildcards:
        if num_unlabeled_wildcards != num_wildcards:
            raise click.UsageError(
                "Internal error. Found mismatching number of unlabeled wildcards "
                f"in --{where} SMILES {smiles!r}")
        return num_wildcards, False
        
    # Otherwise, if there are unlabled wildcards, double-check the numbers match.
    if len(seen_labeled_wildcards) != num_wildcards:
        raise click.UsageError(
            "Internal error. Found mismatching number of labeled wildcards "
            f"in --{where} SMILES {smiles!r}")

    # Check that the correct labels were used
    if "[*:3]" in seen_labeled_wildcards and num_wildcards < 3:
        raise click.UsageError(
            "Cannot have '[*:3]' without specifing both '[*:1]' and '[*:2]' "
            f"in --{where} SMILES {smiles!r}")
        
    if "[*:2]" in seen_labeled_wildcards and num_wildcards < 2:
        raise click.UsageError(
            "Cannot have '[*:2]' without specifing '[*:1]' in "
            f"--{where} SMILES {smiles!r}")
    
    return num_wildcards, True

# Check that the SMILES is valid, with a single fragment, and
# the right wildcard atoms. I've already tokenized the SMILES
# string so I know how many wildcard atoms there should be.
#
# Return the molecule and bonds (as atom pairs) where one end contains
# a wildcard.

wildcard_match_pat = None
def check_valid_smiles(where, smiles, num_expected_wildcards):
    from rdkit import Chem
    
    global wildcard_match_pat
    
    if wildcard_match_pat is None:
        # wildcard with one connection (a non-ring connection) to a non-wildcard, non-hydrogen
        wildcard_match_pat = Chem.MolFromSmarts("[#0X1]!@[!#0,!#1]")

    # Is it a valid SMILES?
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise click.UsageError(f"Cannot parse --{where} SMILES {smiles!r}")

    # Does it contain only one fragment?
    num_frags = len(Chem.GetMolFrags(mol))
    if num_frags != 1:
        raise click.UsageError(
            f"--{where} SMILES {smiles!r} must have one and only one fragment")

    # Look for the wildcard atoms.
    num_frags = len(Chem.GetMolFrags(mol))        
    matches = mol.GetSubstructMatches(wildcard_match_pat)
    if len(matches) != num_expected_wildcards:
        raise click.UsageError(
            f"The wildcard atoms in --{where} SMILES {smiles!r} "
            "do not appear to be valid attachment points")
    return mol, matches


# Check that the "--from" SMILES is correct.

unlabeled_wildcard_pat = re.compile(
    re.escape("[*]") +
    "|" +
    re.escape("*")) # Support both "[*]" and "*" forms

def check_from_smiles(from_smiles):
    from rdkit import Chem
    
    num_wildcards, has_labels = check_wildcards("from", from_smiles)

    # Need to put it into canonical form.
    # TODO: This should be a more readily accessible library function.
    
    if has_labels:
        # Remove the labels.
        from_smiles = (from_smiles
                           .replace("[*:3]", "[*]")
                           .replace("[*:2]", "[*]")
                           .replace("[*:1]", "[*]")
                           )

    # Canonicalize without labels
    mol, cutlist = check_valid_smiles("from", from_smiles, num_wildcards)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    # Assign the values *:1, *:2, and *:3 in order
    # (Listed backwards because I pop() them off.)
    sub_terms = ["[*:3]", "[*:2]", "[*:1]"]
    from_smiles = unlabeled_wildcard_pat.sub(lambda pat: sub_terms.pop(), from_smiles)

    return from_smiles

# The "cansmirks" implementation assumes there are no wildcard atoms.
# I can fake the algorithm by replacing the '*' atoms with an unused element.
def prepare_mol_for_cansmirks(where, mol):
    atomic_nums = set(a.GetAtomicNum() for a in mol.GetAtoms())
    for atomic_num in range(120, 2, -1):
        if atomic_num not in atomic_nums:
            break
    else:
        raise click.UsageError(
            f"Cannot work with --{where} SMILES because it contains every usable element")
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(atomic_num)

# If --from and --to are given, then life is much more complicated.
# Need to also ensure they are compatible, and do different work if
# given labled and unlabeled SMILES.
# Returns a list of unique SMIRKS pairs.
def check_from_to_smiles(from_smiles, to_smiles):
    # First, check that the --from and --to SMILES are reasonable
    num_from_wildcards, from_has_labels = check_wildcards("from", from_smiles)
    num_to_wildcards, to_has_labels = check_wildcards("to", to_smiles)

    # Check they have the same number of wildcards
    if num_from_wildcards != num_to_wildcards:
        raise click.UsageError(
            f"--from SMILES {from_smiles!r} and --to SMILES {to_smiles!r} "
            f"must have the same number of wildcard atoms: "
            f"{num_from_wildcards} != {num_to_wildcards}")

    # Check they are both labeled or both unlabeled
    if from_has_labels != to_has_labels:
        raise click.UsageError(
            "Cannot mix unlabeled atoms in the --from SMILES "
            "with labeled atoms in the --to SMILES")

    if num_from_wildcards == 1:
        # Single cut is easy
        return [(from_smiles, to_smiles)]

    # will cut the wildcard atoms off, and cansmirks the corresponding non-wildcard fragments.
    from_mol, from_cutlist = check_valid_smiles("from", from_smiles, num_from_wildcards)
    prepare_mol_for_cansmirks("from", from_mol)
    to_mol, to_cutlist = check_valid_smiles("to", to_smiles, num_to_wildcards)
    prepare_mol_for_cansmirks("to", to_mol)
        
    if from_has_labels:
        # Both are labeled
        return create_labeled_cansmirks(num_from_wildcards, from_mol, from_cutlist, to_mol, to_cutlist)
    else:
        # Both are unlabeled
        return create_unlabeled_cansmirks(num_from_wildcards, from_mol, from_cutlist, to_mol, to_cutlist)

# Helper class to pass the wildcard cut bonds to fragment_algorithm.fragment_mol()
class FragmentCutlist(object):
    min_heavies_per_const_frag = 0
    min_heavies_total_const_frag = 0
    def __init__(self, cutlist):
        self.cutlist = cutlist
    def get_cut_lists(self, mol):
        return [self.cutlist]


# This is the easy one since the mapping is already specified.
def create_labeled_cansmirks(from_mol, from_cutlist, to_mol, to_cutlist):
    from .. import fragment_algorithm
    from .. import index_algorithm
    
    num_cuts = len(from_cutlist)
    constant_symmetry_class = "123"[:num_cuts]

    # Fragment the --from and --to along the given cutlists
    from_fragmentation = next(
        fragment_algorithm.fragment_mol(from_mol, FragmentCutlist(from_cutlist), 100)
        )
    to_fragmentation = next(
        fragment_algorithm.fragment_mol(to_mol, FragmentCutlist(to_cutlist), 100)
        )
    # Canonicalize the ordering of the two sides
    need_swap = (
        (from_fragmentation.variable_smiles, from_fragmentation.attachment_order) >
        (to_fragmentation.variable_smiles, to_fragmentation.attachment_order)
        )
    if need_swap:
        from_fragmentation, to_fragmentation = to_fragmentation, from_fragmentation
        
    
    # Use 'cansmirks' to figure out the canonical labeling.
    smirks, constant_smiles = index_algorithm.cansmirks(
        num_cuts,
        from_fragmentation.variable_smiles, from_fragmentation.variable_symmetry_class, from_fragmentation.attachment_order,
        "", constant_symmetry_class,
        to_fragmentation.variable_smiles, to_fragmentation.variable_symmetry_class, to_fragmentation.attachment_order,
        index_algorithm.RelabelCache())

    # Return the order to the user-specified one (--from/--to)
    if need_swap:
        return [smirks.split(">>")[::-1]]
    else:
        return [smirks.split(">>")]
        

# The two molecules are unlabeled, so we need to identify up to n! canonical forms.
def create_unlabeled_cansmirks(from_mol, from_cutlist, to_mol, to_cutlist):
    from .. import index_algorithm
    
    num_cuts = len(from_cutlist)
    constant_symmetry_class = "123"[:num_cuts]

    # Fragment the --from and --to along the given cutlists
    from_fragmentation = next(
        fragment_algorithm.fragment_mol(from_mol, FragmentCutlist(from_cutlist), 100)
        )
    to_fragmentation = next(
        fragment_algorithm.fragment_mol(to_mol, FragmentCutlist(to_cutlist), 100)
        )
    
    # Put into canonical order (required for 'cansmirks')
    need_swap = (
        (from_fragmentation.variable_smiles, from_fragmentation.attachment_order) >
        (to_fragmentation.variable_smiles, to_fragmentation.attachment_order)
        )
    if need_swap:
        from_fragmentation, to_fragmentation = to_fragmentation, from_fragmentation

    # Enumerate all n! ways
    if num_cuts == 2:
        from_attachment_orders = ("01", "10")
    elif num_cuts == 3:
        from_attachment_orders = ("012", "021", "120", "102", "210", "201")
    else:
        raise AssertionError(num_cuts) # n=1 is special-cased elsewhere

    unique_smirks = set()
    for from_attachment_order in from_attachment_orders:
        smirks, constant_smiles = index_algorithm.cansmirks(
            num_cuts,
            from_fragmentation.variable_smiles, from_fragmentation.variable_symmetry_class, from_attachment_order,
            "", constant_symmetry_class,
            to_fragmentation.variable_smiles, to_fragmentation.variable_symmetry_class, to_fragmentation.attachment_order,
            index_algorithm.RelabelCache())
        unique_smirks.add(smirks)

    # Return the order to the user-specified one (--from/--to)
    if need_swap:
        return [smirks.split(">>")[::-1] for smirks in sorted(unique_smirks)]
    else:
        return [smirks.split(">>") for smirks in sorted(unique_smirks)]
        

def get_environment_fragments_up_to_radius(mol, max_radius=5):
    assert max_radius >= 0
    atom_indices = set()
    bond_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom_indices.add(atom.GetIdx())

    results = [(0, sorted(atom_indices), sorted(bond_indices))]
            
    old_atom_indices = set(atom_indices)
    for r in range(max_radius):
        new_atom_indices = set()
        for atom_idx in old_atom_indices:
            for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
                bond_idx = bond.GetIdx()
                if bond_idx in bond_indices:
                    continue
                bond_indices.add(bond_idx)
                other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
                if other_atom_idx in atom_indices:
                    continue
                new_atom_indices.add(other_atom_idx)
        if not new_atom_indices:
            break
        atom_indices.update(new_atom_indices)
        old_atom_indices = new_atom_indices

        results.append( (r+1, sorted(atom_indices), sorted(bond_indices)) )
        
    return results

def get_environment_smiles_up_to_radius(smiles, max_radius=5):
    from rdkit import Chem
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [(r, "<bad smiles>") for r in range(max_radius+1)]
    results = []
    for (r, atom_indices, bond_indices) in get_environment_fragments_up_to_radius(mol, max_radius):
        fragment_smiles = Chem.MolFragmentToSmiles(mol, atom_indices, bond_indices)
        results.append( (r, fragment_smiles) )
    return results

class EnvironmentSmilesLookup(object):
    def __init__(self, dataset):
        self.dataset = dataset
        self._cache = {}
        
    def get_environment_smiles(self, rule_environment_id, radius, cursor):
        key = (rule_environment_id, radius)
        if key in self._cache:
            return self._cache[key]
        
        constant_smiles = self.dataset.get_smallest_constant_smiles(rule_environment_id, cursor)

        for (r, fragment_smiles) in get_environment_smiles_up_to_radius(constant_smiles):
            self._cache[(rule_environment_id, r)] = fragment_smiles
            
        return self._cache.get(key, "<missing>") # XXX: shouldn't this always be present?
        
    
####

@command()

@click.option(
    "--from",
    "from_smiles",
    metavar = "SMILES",
    default = None,
    help = "SMILES for one side of the transformation",
    )

@click.option(
    "--to",
    "to_smiles",
    metavar = "SMILES",
    default = None,
    help = "SMILES for the other side of the transformation",
    )

@click.option(
    "--canonicalize/--no-canonicalize",
    default = True,
    help = "Use the --from and --to strings as-is; do not canonicalize them (default: --canonicalize)",
    )

@click.option(
    "--property",
    "-p",
    "property_names",
    metavar="NAME",
    multiple=True,
    help="Property to use (may be specified multiple times)",
    )

@click.option(
    "--min-count",
    type = positive_int(),
    help = "Only show rules with at least N pairs",
    )

@click.option(
    "--output",
    "-o",
    "outfile",
    default = "-",
    type = GzipFile("w"),
    help = "Write the output to the given file (default is stdout)",
    )


@add_single_database_parameters()

@click.pass_obj
def proprulecat(
        reporter,
        from_smiles,
        to_smiles,
        canonicalize,
        property_names,
        min_count,
        outfile,
        database_options,
        ):
    """Write the property rules to stdout or a file"""

    dataset = open_dataset_from_options_or_exit(database_options)
    c = dataset.get_cursor()

    # Validate the property names and get the corresponding property ids
    property_names = get_property_names_or_error(
        dataset,
        property_names = property_names,
        all_properties = True,
        )
    
    property_searches = []
    for property_name in property_names:
        property_searches.append( (dataset.get_property_name_id(property_name), property_name) )
    
    assert min_count is None or min_count > 0, min_count

    if from_smiles is None and to_smiles is not None:
        raise click.UsageError("--from SMILES is required to be able to specify a --to SMILES")

    smirks_pairs = []

    if not canonicalize:
        smirks_pairs = [(from_smiles, to_smiles)]
        
    elif from_smiles is None and to_smiles is None:
        smirks_pairs = [(None, None)]
        
    elif to_smiles is None:
        # Only need to handle the --from SMILES
        from_smiles = check_from_smiles(from_smiles)
        smirks_pairs.append( (from_smiles, None) )
        
    else:
        # Handle both --from and --to SMILES
        smirks_pairs.extend(check_from_to_smiles(from_smiles, to_smiles))

    
    # If the field is None, then replace it with the "*"
    def NULLABLE(x):
        if x is None:
            return "*"
        return x

    fragment_cursor = dataset.get_cursor()
    environment_smiles_lookup = EnvironmentSmilesLookup(dataset)
    
    with outfile:
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                          ("from_smiles", "to_smiles", "property", "radius", "environment",
                               "count", "avg", "std", "kurtosis", "skewness",
                          "min", "q1", "median", "q3", "max", "paired_t", "p_value")
                          )
    
        for property_id, property_name in property_searches:
            for from_smiles, to_smiles in smirks_pairs:
                property_rules_iter = dataset.iter_selected_property_rules(
                    from_smiles, to_smiles, property_id, min_count=min_count, cursor=c)
                
                for property_rule in property_rules_iter:
                    R = property_rule

                    fragment_smiles = environment_smiles_lookup.get_environment_smiles(
                        R.rule_environment_id, R.radius, fragment_cursor)
                    
                    if min_count is not None and R.count < min_count:
                        raise AssertionError("min_count failed in the database search")
                    outfile.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        R.from_smiles, R.to_smiles, property_name, R.radius,
                        fragment_smiles,
                        R.count, R.avg, NULLABLE(R.std), NULLABLE(R.kurtosis), NULLABLE(R.skewness),
                        R.min, R.q1, R.median, R.q3, R.max, NULLABLE(R.paired_t), NULLABLE(R.p_value)
                    ))

    
