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
import gzip
import re

from mmpdblib import commandline

from support import capture_stdout, capture_stderr
from support import get_filename, create_test_filename

from rdkit import Chem
wildcard_atom = Chem.CanonSmiles("*")
if wildcard_atom == "[*]":
    TEST_DATA_MMPDB = get_filename("test_data_2019.mmpdb")
elif wildcard_atom == "*":
    # The dataset was generated with:
    #   python -m mmpdblib.commandline index test_data.fragments --properties test_data.csv -o test_data_2018.mmpdb
    TEST_DATA_MMPDB = get_filename("test_data_2019.mmpdb")
else:
    raise AssertionError(wildcard_atom)
    

class Table(list):
    def get_column(self, column_name):
        i = self[0].index(column_name)
        return [row[i] for row in self[1:]]

def parse_table(output):
    rows = [line.split("\t") for line in output.splitlines(False)]
    return Table(rows)

def transform(*args):
    args = ("--quiet", "transform", TEST_DATA_MMPDB) + tuple(args)
    with capture_stdout() as stdout:
        try:
            commandline.main(args)
        except SystemExit as err:
            raise AssertionError("SystemExit trying to run %r: %s" % (args, err))
    return parse_table(stdout.value)

def transform_fail(*args):
    args = ("--quiet", "transform", TEST_DATA_MMPDB) + tuple(args)
    with capture_stderr() as stderr:
        try:
            commandline.main(args)
        except SystemExit as err:
            pass
        else:
            raise AssertionError("Should have failed: %r" % (args,))
    return stderr.value
    

class TestTransformCommand(unittest.TestCase):
    def check_property_delta(self, table, property, expected_dict):
        expected_dict = expected_dict.copy()
        for from_smiles, to_smiles, delta in zip(
                table.get_column(property + "_from_smiles"),
                table.get_column(property + "_to_smiles"),
                table.get_column(property + "_avg")):
            self.assertEqual(from_smiles, "[*:1]O")
            expected = expected_dict.pop(to_smiles)
            self.assertEqual(float(delta), expected)
        self.assertEqual(expected_dict, {})

    def test_basic(self):
        table = transform("--smiles", "c1cccnc1O")
        self.assertEqual(table[0], [
            'ID',
            'SMILES',
            'MW_from_smiles',
            'MW_to_smiles',
            'MW_radius',
            'MW_fingerprint',
            'MW_rule_environment_id',
            'MW_count',
            'MW_avg',
            'MW_std',
            'MW_kurtosis',
            'MW_skewness',
            'MW_min',
            'MW_q1',
            'MW_median',
            'MW_q3',
            'MW_max',
            'MW_paired_t',
            'MW_p_value',
            'MP_from_smiles',
            'MP_to_smiles',
            'MP_radius',
            'MP_fingerprint',
            'MP_rule_environment_id',
            'MP_count',
            'MP_avg',
            'MP_std',
            'MP_kurtosis',
            'MP_skewness',
            'MP_min',
            'MP_q1',
            'MP_median',
            'MP_q3',
            'MP_max',
            'MP_paired_t',
            'MP_p_value',
            ])

        # Make sure the changes are as expected
        self.check_property_delta(table, "MW", {
            # remember, *O means *[OH] but *Cl is just *[Cl]
            '[*:1]Cl': 35.5 - 16.0 - 1.0,
            # *N is *[NH2]
            '[*:1]N': 14.0 + 2.0 - 16.0 - 1.0,
            # [*:1][H] is just a hydrogen
            '[*:1][H]': 1.0 - 16.0 - 1.0})

        self.check_property_delta(table, "MP", {
            '[*:1]Cl': -97.0,
            '[*:1]N': -16.667,
            '[*:1][H]': -93.0})

        self.assertEqual(table.get_column("MW_count"), ["1", "3", "4"])
        self.assertEqual(table.get_column("MP_count"), ["1", "3", "3"])

        self.assertEqual(table.get_column("MW_radius"), ["1", "1", "1"])

    def test_MP_property(self):
        table = transform("--smiles", "c1cccnc1O", "--property", "MP")
        self.assertEqual(table[0], [
            'ID',
            'SMILES',
            'MP_from_smiles',
            'MP_to_smiles',
            'MP_radius',
            'MP_fingerprint',
            'MP_rule_environment_id',
            'MP_count',
            'MP_avg',
            'MP_std',
            'MP_kurtosis',
            'MP_skewness',
            'MP_min',
            'MP_q1',
            'MP_median',
            'MP_q3',
            'MP_max',
            'MP_paired_t',
            'MP_p_value',
            ])

        # Make sure the changes are as expected
        self.check_property_delta(table, "MP", {
            '[*:1]Cl': -97.0,
            '[*:1]N': -16.667,
            '[*:1][H]': -93.0})

    def test_no_properties(self):
        table = transform("--smiles", "c1cccnc1O", "--no-properties")
        self.assertEqual(table[0], ["ID", "SMILES"])
        self.assertEqual(len(table), 4)

    def test_radius(self):
        # radius=1 doesn't change things because that's still the aromatic carbons
        table = transform("--smiles", "c1cccnc1O", "--min-radius", "1", "--property", "MW")
        self.check_property_delta(table, "MW", {
            '[*:1]Cl': 35.5 - 16.0 - 1.0,
            '[*:1]N': 14.0 + 2.0 - 16.0 - 1.0,
            '[*:1][H]': 1.0 - 16.0 - 1.0})

        # radius=2 reaches the aromatic nitrogen, so finds no matches
        table = transform("--smiles", "c1cccnc1O", "--min-radius", "2", "--property", "MW")
        self.check_property_delta(table, "MW", {})

    def test_substructure(self):
        table = transform("--smiles", "c1cccnc1O", "--substructure", "Cl", "--property", "MW")
        self.check_property_delta(table, "MW", {
            '[*:1]Cl': 35.5 - 16.0 - 1.0})

    def test_min_pairs(self):
        table = transform("--smiles", "c1cccnc1O", "--min-pairs", "3")
        self.assertEqual(table.get_column("MW_count"), ["3", "4"])
        self.assertEqual(table.get_column("MP_count"), ["3", "3"])

        table = transform("--smiles", "c1cccnc1O", "--min-pairs", "4")
        self.assertEqual(table.get_column("MW_count"), ["4"])
        self.assertEqual(table.get_column("MP_count"), [""])

    def test_where(self):
        # This is the more complicated way to specify a --min-pairs filter.
        table = transform("--smiles", "c1cccnc1O", "--where", "count > 2")
        self.assertEqual(table.get_column("MW_count"), ["3", "4"])
        self.assertEqual(table.get_column("MP_count"), ["3", "3"])

        table = transform("--smiles", "c1cccnc1O", "--where", "count > 3")
        self.assertEqual(table.get_column("MW_count"), ["4"])
        self.assertEqual(table.get_column("MP_count"), [""])

    def test_score(self):
        # Choose the smallest radius (note the leading space so it doesn't get
        # confused for a command-line option)
        table = transform("--smiles", "c1cccnc1O", "--score", " -min-radius", "--property", "MW")
        self.assertEqual(table.get_column("MW_radius"), ["0", "0", "0"])

    def test_radius(self):
        # Require a specific radius
        table = transform("--smiles", "c1cccnc1O", "--min-radius", "2")
        self.assertEqual(table.get_column("MW_radius"), [])
        table = transform("--smiles", "c1cccnc1O", "--min-radius", "1")
        self.assertEqual(table.get_column("MW_radius"), ["1", "1", "1"])

    # XXX I have no test for --rule-selection-cutoffs. Does my test set support it?

    def test_output(self):
        filename = create_test_filename(self, "transform.out")
        transform("--smiles", "c1cccnc1O", "--output", filename)
        with open(filename) as infile:
            table = parse_table(infile.read())
        self.assertEqual(len(table), 4)
        self.assertIn("MP_count", table[0])
        self.assertIn("MW_count", table[0])

    def test_output_gz(self):
        filename = create_test_filename(self, "transform.out.gz")
        transform("--smiles", "c1cccnc1O", "--output", filename)
        with gzip.open(filename) as infile:
            table = parse_table(infile.read().decode("ascii"))
        self.assertEqual(len(table), 4)
        self.assertIn("MP_count", table[0])
        self.assertIn("MW_count", table[0])

    #### Error conditions
    def test_bad_smiles(self):
        errmsg = transform_fail("--smiles", "NOT_A_SMILES")
        self.assertIn("error: Unable to fragment --smiles 'NOT_A_SMILES': invalid smiles\n", errmsg)

    def test_bad_smiles_not_enough_heavy_atoms(self):
        errmsg = transform_fail("--smiles", "C")
        self.assertIn("error: Unable to fragment --smiles 'C': not enough heavy atoms\n", errmsg)

    def test_bad_min_constant_size(self):
        #[--min-constant-size N]
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--min-constant-size", "Q")
        self.assertIn("--min-constant-size: must be a positive integer or zero\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--min-constant-size", "-1")
        self.assertIn("--min-constant-size: must be a positive integer or zero\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--min-constant-size", "1.0")
        self.assertIn("--min-constant-size: must be a positive integer or zero\n", errmsg)

    def test_bad_radius(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--min-radius", "-1")
        self.assertIn("--min-radius/-r: invalid choice: '-1'", errmsg)
        self.assertIn("'0', '1', '2', '3', '4', '5'", errmsg)

    def test_bad_min_pairs(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--min-pairs", "-1")
        self.assertIn("--min-pairs: must be a positive integer or zero\n", errmsg)

    def test_bad_substructure(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--substructure", "ZZTop")
        self.assertIn("Cannot parse --substructure 'ZZTop'\n", errmsg)

    def test_bad_property_name(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--property", "BP")
        self.assertIn("--property 'BP' is not present in the database\n", errmsg)

    def test_bad_where(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--where", "BAD_VARIABLE")
        self.assertIn("unsupported variable name 'BAD_VARIABLE', in --where 'BAD_VARIABLE'\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--where", "invalid Python expression")
        self.assertIn("invalid syntax (--where, line 1), in --where 'invalid Python expression'\n", errmsg)

    def test_bad_score(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--score", "BAD_VARIABLE")
        self.assertIn("unsupported variable name 'BAD_VARIABLE', in --score 'BAD_VARIABLE'\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--score", "invalid Python expression")
        self.assertIn("invalid syntax (--score, line 1), in --score 'invalid Python expression'\n", errmsg)

    def test_bad_rule_selection_cutoffs(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--rule-selection-cutoffs", "A,B,C")
        self.assertIn("--rule-selection-cutoffs: could not parse 'A' as an integer: "
                      "invalid literal for int() with base 10: 'A'\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--rule-selection-cutoffs", "20,10,-4")
        self.assertIn("--rule-selection-cutoffs: threshold values must be non-negative\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--rule-selection-cutoffs", "10,20")
        self.assertIn("--rule-selection-cutoffs: threshold values must be in decreasing order\n", errmsg)

    def test_bad_jobs(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--jobs", "-1")
        self.assertIn("--jobs/-j: must be a positive integer\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--jobs", "0")
        self.assertIn("--jobs/-j: must be a positive integer\n", errmsg)
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--jobs", "0 0")
        self.assertIn("--jobs/-j: must be a positive integer\n", errmsg)

    def test_bad_output(self):
        errmsg = transform_fail("--smiles", "c1cccnc1O", "--output", "/this/directory/does/not/at/all/exist/output.txt")
        self.assertIn("Cannot open --output file:", errmsg)
        self.assertIn("/this/directory/does/not/at/all/exist/output.txt", errmsg)

######

_predict_pat = re.compile(r"^predicted delta: (\S+)( predicted value: (\S+))?$")

def parse_predict(text):
    if text == "No rules found.\n":
        return None, None
        
    m = _predict_pat.match(text)
    if m is None:
        raise AssertionError("Could not parse %r" % (text,))
    return m.group(1), m.group(3)
        
def predict(*args):
    args = ("--quiet", "predict", TEST_DATA_MMPDB) + tuple(args)
    with capture_stdout() as stdout:
        try:
            commandline.main(args)
        except SystemExit as err:
            raise AssertionError("SystemExit trying to run %r: %s" % (args, err))
    return parse_predict(stdout.value)

def predict_fail(*args):
    args = ("--quiet", "predict", TEST_DATA_MMPDB) + tuple(args)
    with capture_stderr() as stderr:
        try:
            commandline.main(args)
        except SystemExit as err:
            pass
        else:
            raise AssertionError("Should have failed: %r" % (args,))
    return stderr.value

def read_details(prefix):
    with open(prefix + "_rules.txt") as infile:
        rules_table = parse_table(infile.read())
    with open(prefix + "_pairs.txt") as infile:
        pairs_table = parse_table(infile.read())
    return rules_table, pairs_table

class TestPredictCommand(unittest.TestCase):
    def test_MW_basic(self):
        delta, new_value = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                                   "--property", "MW")
        self.assertEqual(delta, "-18.5")
        self.assertIs(new_value, None)

    def test_MW_value(self):
        delta, new_value = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                                   "--property", "MW", "--value", "20.1")
        self.assertEqual(delta, "-18.5")
        self.assertEqual(new_value, "1.6")

    def test_MP_basic(self):
        delta, new_value = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                                   "--property", "MP")
        self.assertEqual(delta, "+97")
        self.assertIs(new_value, None)

    def test_MP_value(self):
        delta, new_value = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                                   "--property", "MP", "--value", "1.23")
        self.assertEqual(delta, "+97")
        self.assertEqual(new_value, "98.23")

    def test_where(self):
        result = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                         "--property", "MP", "--where", "count > 10")
        self.assertEqual(result, (None, None))

    def test_save_details(self):
        prefix = create_test_filename(self, "predict_details")
        result = predict("--smiles", "c1cccnc1O", "--reference", "Clc1ccccn1",
                         "--property", "MW", "--save-details",
                         "--prefix", prefix)
        self.assertEqual(result, ("-18.5", None))

        rules_table, pairs_table = read_details(prefix)
        self.assertEqual(rules_table.get_column("from_smiles"), ["[*:1]Cl", "[*:1]Cl"])
        self.assertEqual(rules_table.get_column("to_smiles"), ["[*:1]O", "[*:1]O"])
        
        self.assertEqual(pairs_table.get_column("radius"), ["0", "1"])
        self.assertEqual(pairs_table.get_column("lhs_public_id"), ["2-chlorophenol", "2-chlorophenol"])
        self.assertEqual(pairs_table.get_column("rhs_public_id"), ["catechol", "catechol"])
        

    
        
if __name__ == "__main__":
    unittest.main()
