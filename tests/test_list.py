import unittest
import re
import support
from support import get_filename

TEST_DATA_MMPDB = get_filename("test_data_2019.mmpdb")

def list_main(args):
    args = ("--quiet", "list") + tuple(args)
    return support.expect_pass(args)


def list_main_fail(args):
    args = ("--quiet", "list") + tuple(args)
    return support.expect_fail(args).stderr


_header_fields = """
Name #cmpds #rules #pairs #envs  #stats |Title| Properties
""".split()

_expected_fields = f"""
{TEST_DATA_MMPDB}      9     47    342    321    533  MMPs from 'test_data.fragdb' MW MP
""".split()


class TestList(unittest.TestCase):
    def _check_header(self, line):
        fields = []
        # Normalize the '|--- Title ---|' section
        line = re.sub(r"\|-{1,} ", "|", line)
        line = re.sub(r" -{1,}\|", "|", line)
        for field in line.split():
            fields.append(field)

        self.assertEqual(fields, _header_fields)

    def _check_output(self, output):
        lines = output.splitlines()
        self._check_header(lines[0])
        n = 1
        if len(lines) > 1:
            n = 0 
            for line in lines[1:]:
                if "test_data_2019.mmpdb" in line:
                    n += 1
                    self.assertEqual(line.split(), _expected_fields)
        return n

    def test_no_args(self):
        result = list_main([])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_one_arg(self):
        result = list_main([TEST_DATA_MMPDB])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_two_args(self):
        result = list_main([TEST_DATA_MMPDB, TEST_DATA_MMPDB])
        n = self._check_output(result.output)
        self.assertEqual(n, 2, "should have existed twice")

    def test_file_does_not_exist(self):
        result = list_main(["does_not_exist.mmpdb", TEST_DATA_MMPDB])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_recount(self):
        result = list_main(["--recount", TEST_DATA_MMPDB])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_all(self):
        self._test_all(["--all", TEST_DATA_MMPDB])

    def test_a(self):
        self._test_all(["-a", TEST_DATA_MMPDB])

    def _test_all(self, flag):
        result = list_main(flag)
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

        # If the data set is ever regenerated then at the very least
        # the 'Created' line will need to be updated.
        expected_lines = f"""\
                         Name                         #cmpds #rules #pairs #envs  #stats  |--------- Title ----------| Properties
{TEST_DATA_MMPDB}      9     47    342    321    533  MMPs from 'test_data.fragdb' MW MP
      Created: 2025-05-02 14:54:33.639458
        #compounds/property:  8/MP 9/MW
        #smiles for rules: 21  for constants: 10
        Fragment options:
          cut_smarts: [#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]
          max_heavies: 100
          max_rotatable_bonds: 10
          max_up_enumerations: 1000
          method: chiral
          min_heavies_per_const_frag: 0
          min_heavies_total_const_frag: 0
          num_cuts: 3
          rotatable_smarts: [!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]
          salt_remover: <default>
        Index options:
          max_radius: 5
          max_variable_heavies: 10
          min_radius: 0
          smallest_transformation_only: False
          symmetric: False""".splitlines()
        result_lines = result.output.splitlines()
        num_checked = 0
        if len(expected_lines) == len(result_lines):
            for line1, line2 in zip(expected_lines, result_lines):
                self.assertIn(re.sub(r'\s+', ' ', line1).strip(), re.sub(r'\s+', ' ', line2).strip())
                num_checked += 1
        self.assertEqual(num_checked, 22)


if __name__ == "__main__":
    unittest.main()
