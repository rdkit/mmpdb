import unittest
import support


def list_main(args):
    args = ("--quiet", "list") + tuple(args)
    return support.expect_pass(args)


def list_main_fail(args):
    args = ("--quiet", "list") + tuple(args)
    return support.expect_fail(args).stderr


_header_fields = """
Name         #cmpds #rules #pairs #envs  #stats  |- Title -| Properties
""".split()

_expected_fields = """\
test_data_2019.mmpdb      9     47    342    321    533  MMPs from 'test_data.fragdb'                                    MW MP
""".split()


class TestList(unittest.TestCase):
    def _check_header(self, line):
        fields = []
        for field in line.split():
            # Normalize the '|--- Title ---|' section
            if field.startswith("|-"):
                field = "|-"
            elif field.endswith("-|"):
                field = "-|"
            fields.append(field)

        self.assertEqual(fields, _header_fields)

    def _check_output(self, output):
        lines = output.splitlines()
        self._check_header(lines[0])
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
        result = list_main(["test_data_2019.mmpdb"])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_two_args(self):
        result = list_main(["test_data_2019.mmpdb", "test_data_2019.mmpdb"])
        n = self._check_output(result.output)
        self.assertEqual(n, 2, "should have existed twice")

    def test_file_does_not_exist(self):
        result = list_main(["does_not_exist.mmpdb", "test_data_2019.mmpdb"])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_recount(self):
        result = list_main(["--recount"])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

    def test_all(self):
        self._test_all("--all")

    def test_a(self):
        self._test_all("-a")

    def _test_all(self, flag):
        result = list_main([flag])
        n = self._check_output(result.output)
        self.assertEqual(n, 1)

        # If the data set is ever regenerated then at the very least
        # the 'Created' line will need to be updated.
        num_checked = 0
        for expected_line in """
        Name         #cmpds #rules #pairs #envs  #stats  |--------- Title ----------| Properties
test_data_2019.mmpdb      9     47    342    321    533  MMPs from 'test_data.fragdb' MW MP
      Created: 2021-12-03 13:57:35.268963
        #compounds/property:  8/MP 9/MW
        #smiles for rules: 21  for constants: 10
        Fragment options:
          cut_smarts: [#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]
          max_heavies: 100
          max_rotatable_bonds: 10
          method: chiral
          min_heavies_per_const_frag: 0
          num_cuts: 3
          rotatable_smarts: [!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]
          salt_remover: <default>
        Index options:
          max_radius: 5
          max_variable_heavies: 10
          min_radius: 0
          smallest_transformation_only: False
          symmetric: False
        """.splitlines():
            expected_line = expected_line.strip()
            if not expected_line:
                continue
            self.assertIn(expected_line, result.output)
            num_checked += 1

        self.assertEqual(num_checked, 20)


if __name__ == "__main__":
    unittest.main()
