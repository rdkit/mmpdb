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

import unittest
import json
import sqlite3
from contextlib import closing

from mmpdblib import config, dbutils
from mmpdblib.index_algorithm import (
    EnvironmentCache,
    VariableFragmentsReducer,
    get_max_constant_radius_for_fraction_transfer,
)

from support import get_filename, create_test_filename, expect_pass

TEST_DATA_FRAGDB = get_filename("test_data.fragdb")
TEST_DATA_CSV = get_filename("test_data.csv")


def index(mmpdb_filename, *args):
    args = ("--quiet", "index", TEST_DATA_FRAGDB, "-o", mmpdb_filename) + tuple(args)
    expect_pass(args)


class TestIndexCommandline(unittest.TestCase):
    def _get_options(self, *args):
        mmpdb_filename = create_test_filename(self, "default.mmpdb")
        index(mmpdb_filename, *args)
        db = dbutils.open_database(mmpdb_filename)
        dataset = db.get_dataset()
        return dataset, json.loads(dataset.index_options_str)

    def test_default_filters(self):
        dataset, options = self._get_options()
        self.assertEqual(
            options,
            {
                "min_radius": 0,
                "max_radius": 5,
                "max_variable_heavies": 10,
                "smallest_transformation_only": False,
                "symmetric": False,
            },
        )
        self.assertEqual(dataset.get_num_rules(), 47)
        self.assertEqual(dataset.get_num_pairs(), 342)
        self.assertEqual(dataset.get_property_names(), [])

    def test_set_filters(self):
        dataset, options = self._get_options(
            "--min-variable-heavies",
            "1",
            "--max-variable-heavies",
            "29",
            "--min-variable-ratio",
            "0.1",
            "--max-variable-ratio",
            "0.99",
            "--max-heavies-transf",
            "25",
            "--symmetric",
            "--max-frac-trans",
            "3",
        )
        self.assertEqual(
            options,
            {
                "symmetric": True,
                "max_frac_trans": 3.0,
                "max_heavies_transf": 25,
                "min_radius": 0,
                "max_radius": 5,
                "max_variable_heavies": 29,
                "max_variable_ratio": 0.99,
                "min_variable_heavies": 1,
                "min_variable_ratio": 0.1,
                "smallest_transformation_only": False,
            },
        )
        self.assertEqual(dataset.get_num_rules(), 2 * 47)  # because --symmetric
        self.assertEqual(dataset.get_num_pairs(), 2 * 342)
        self.assertEqual(dataset.get_property_names(), [])

    def test_max_variable_heavies_none(self):
        dataset, options = self._get_options("--max-variable-heavies", "none")
        self.assertEqual(
            options,
            {
                "min_radius": 0,
                "max_radius": 5,
                "symmetric": False,
                "smallest_transformation_only": False,
            },
        )

    def test_with_properties(self):
        dataset, options = self._get_options("--properties", TEST_DATA_CSV, "--title", "test data")
        self.assertEqual(dataset.title, "test data")
        self.assertEqual(dataset.get_property_names(), ["MW", "MP"])

    def test_smallest_transformation_only(self):
        dataset, options = self._get_options("--smallest-transformation-only")
        self.assertTrue(options["smallest_transformation_only"])
        self.assertEqual(dataset.get_num_rules(), 36)
        self.assertEqual(dataset.get_num_pairs(), 252)


def get_pair_radii(mmpdb_filename, min_radius=0):
    # Return the (compound1, compound2, radius) for each pair environment with radius >= min_radius
    with closing(sqlite3.connect(mmpdb_filename)) as conn:
        return set(
            conn.execute(
                """
                SELECT c1.public_id, c2.public_id, rule_environment.radius
                  FROM pair, compound AS c1, compound AS c2, rule_environment
                 WHERE pair.compound1_id = c1.id
                   AND pair.compound2_id = c2.id
                   AND pair.rule_environment_id = rule_environment.id
                   AND rule_environment.radius >= ?
                """,
                (min_radius,),
            )
        )


class TestIndexMinRadiusWithMaxFracTrans(unittest.TestCase):
    # A best radius below --min-radius used to give rules without any
    # environment, or an AssertionError when computing the environments.
    def _index(self, *args):
        mmpdb_filename = create_test_filename(self, "min_radius.mmpdb")
        index(mmpdb_filename, *args)
        return mmpdb_filename

    def _get_num_rules_and_pairs(self, mmpdb_filename):
        db = dbutils.open_database(mmpdb_filename)
        try:
            dataset = db.get_dataset()
            return dataset.get_num_rules(), dataset.get_num_pairs()
        finally:
            db.close()

    def test_no_rules_without_environments(self):
        mmpdb_filename = self._index("--min-radius", "2", "--max-frac-trans", "0.3")
        self.assertEqual(self._get_num_rules_and_pairs(mmpdb_filename), (0, 0))

    def test_best_radius_below_min_radius(self):
        mmpdb_filename = self._index("--min-radius", "3", "--max-frac-trans", "0.5")
        self.assertEqual(self._get_num_rules_and_pairs(mmpdb_filename), (0, 0))

    def test_same_environments_as_min_radius_0(self):
        expected = get_pair_radii(self._index("--max-frac-trans", "0.5"), min_radius=1)
        self.assertEqual(len(expected), 23)
        pair_radii = get_pair_radii(self._index("--min-radius", "1", "--max-frac-trans", "0.5"))
        self.assertEqual(pair_radii, expected)


class TestMaxConstantRadiusForFractionTransfer(unittest.TestCase):
    # The constant has 11 heavies, the variable has 1, so N = 12.
    # The number of heavies within radius 0..5 is 0, 1, 3, 5, 6, 7.
    CONSTANT_SMILES = "[*:1]c1ccc(CCN(C)C)cc1"
    SMIRKS = "[*:1]F>>[*:1]Cl"

    def _get_radius(self, max_frac_trans, min_radius, max_radius=5):
        return get_max_constant_radius_for_fraction_transfer(
            max_frac_trans, self.SMIRKS, self.CONSTANT_SMILES, min_radius, max_radius, EnvironmentCache()
        )

    def test_min_radius_0(self):
        self.assertEqual(self._get_radius(0.2, 0), 1)  # (1+1)/12 <= 0.2 < (1+3)/12
        self.assertEqual(self._get_radius(0.5, 0), 3)  # (1+5)/12 <= 0.5 < (1+6)/12
        self.assertEqual(self._get_radius(0.9, 0), 5)

    def test_min_radius_not_reached(self):
        self.assertIsNone(self._get_radius(0.2, 2))
        self.assertIsNone(self._get_radius(0.5, 4))

    def test_min_radius_reached(self):
        # The radius must not depend on min_radius, as long as it is reachable
        for min_radius in range(4):
            self.assertEqual(self._get_radius(0.5, min_radius), 3, min_radius)
        self.assertEqual(self._get_radius(0.9, 5), 5)

    def test_max_frac_trans_1(self):
        self.assertEqual(self._get_radius(1.0, 0), 5)
        self.assertEqual(self._get_radius(1.0, 2, max_radius=4), 4)


class TestVariableFragmentsReducer(unittest.TestCase):
    def setUp(self):
        self.reducer = VariableFragmentsReducer(config.DEFAULT_FRAGMENT_OPTIONS.get_fragment_filter())

    def assertReducible(self, smirks):
        self.assertTrue(self.reducer.transformation_is_reducible(smirks), smirks)

    def assertNotReducible(self, smirks):
        self.assertFalse(self.reducer.transformation_is_reducible(smirks), smirks)

    def test_reducible(self):
        self.assertReducible("[*:1]c1ccccc1>>[*:1]c1ccc(F)cc1")  # H>>F on the phenyl
        self.assertReducible("[*:1]c1ccc(Cl)cc1>>[*:1]c1ccc(F)cc1")  # Cl>>F on the phenyl
        self.assertReducible("[*:1]CCOC>>[*:1]CCOCC")  # C>>CC on the CCO

    def test_not_reducible(self):
        self.assertNotReducible("[*:1]C(=O)O>>[*:1]c1nn[nH]n1")
        self.assertNotReducible("[*:1]Cl>>[*:1]F")

    def test_multiple_cuts_on_one_atom_are_kept(self):
        # Reducing these would lose the chirality of the atom with the cuts
        self.assertNotReducible("[*:1]C([*:2])c1ccccc1>>[*:1]C([*:2])c1ccncc1")
        self.assertNotReducible("[*:1]C([*:2])c1ccccc1>>[*:1]C([*:2])C1CC1")
        self.assertNotReducible("[*:1]C([*:2])([*:3])c1ccccc1>>[*:1]C([*:2])([*:3])c1ccncc1")


if __name__ == "__main__":
    unittest.main()
