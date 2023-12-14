import unittest

from rdkit import Chem

from mmpdblib import fragment_algorithm
from test_fragment import FIX


class TestFragmentAlgorithm(unittest.TestCase):

    def test_fragment_molecule_on_explicit_hydrogens(self):
        smiles = '[H][CH2]O[H]'
        frags = fragment_algorithm.fragment_molecule_on_explicit_hydrogens(smiles)
        keys = [FIX(frag.get_unique_key()) for frag in frags]
        self.assertEqual(keys, ['0.*[H].*CO', '0.*[H].*OC'])


if __name__ == "__main__":
    unittest.main()
