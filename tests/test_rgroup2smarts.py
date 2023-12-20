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

from __future__ import print_function

import sys
import unittest

from mmpdblib import rgroup2smarts

import support


def run(cmd, stderr_ok=False, input=None):
    if isinstance(cmd, str):
        cmd = cmd.split()
    result = support.expect_pass(cmd, input=input)

    if not stderr_ok and result.stderr:
        print(stderr.value, file=sys.stderr)
        raise AssertionError(("unexpected stderr", cmd, stderr.value))
    return result.output, result.stderr


def run_fail(cmd):
    if isinstance(cmd, str):
        cmd = cmd.split()

    return support.expect_fail(cmd).stderr


class TestSmilesOnCommandline(unittest.TestCase):
    def test_one_smiles(self):
        stdout, stderr = run("rgroup2smarts --cut-rgroup *c1ccccc1O")
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2])]\n")

    def test_two_smiles(self):
        stdout, stderr = run("rgroup2smarts --cut-rgroup *c1ccccc1O --cut-rgroup *F")
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]\n")

    def test_four_smiles(self):
        stdout, stderr = run(
            "rgroup2smarts --cut-rgroup *c1ccccc1O --cut-rgroup *F " "--cut-rgroup *Cl --cut-rgroup *[OH]"
        )
        self.assertEqual(
            stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1]),$([ClH0v1]),$([OHv2])]\n"
        )

    def test_smiles_with_isotope_and_charge(self):
        stdout, stderr = run("rgroup2smarts --cut-rgroup *c1ccccc1[16O] --cut-rgroup *-[C+](=O)[O-]")
        self.assertEqual(
            stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]),$([C+H0v4](=[OH0v2])-[O-H0v1])]\n"
        )


def fix_stderr(test_case, stderr, filename):
    first_line, mid, rest = stderr.partition("\n")
    test_case.assertIn(repr(filename), first_line)
    first_line = first_line.replace(repr(filename), "frags.smi")

    return first_line + mid + rest


merged_test_cases = (
    (["*c1ccccc1O"], "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2])]"),
    (["*c1ccccc1O", "*F"], "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]"),
    (
        ["*c1ccccc1O", "*F", "*Cl", "*[OH]"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1]),$([ClH0v1]),$([OHv2])]",
    ),
    (
        ["*c1ccccc1[16O]", "*-[C+](=O)[O-]"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]),$([C+H0v4](=[OH0v2])-[O-H0v1])]",
    ),
)

simple_test_cases = {
    "*c1ccccc1O": "*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]",
    "*F": "*-!@[FH0v1]",
    "*Cl": "*-!@[ClH0v1]",
    "*[OH]": "*-!@[OHv2]",
    "*c1ccccc1[16O]": "*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]",
    "*-[C+](=O)[O-]": "*-!@[C+H0v4](=[OH0v2])-[O-H0v1]",
}


class TestSmilesOnCommandline(unittest.TestCase):
    def _test_merged(self, check):
        for inputs, output in merged_test_cases:
            args = ["rgroup2smarts"]
            if check:
                args.append("--check")
            for smiles in inputs:
                args.extend(["--cut-rgroup", smiles])
            stdout, stderr = run(args)
            self.assertEqual(stdout, output + "\n", repr(inputs))

    def test_merged(self):
        self._test_merged(False)

    def test_merged_check(self):
        self._test_merged(True)

    def _test_single(self, check):
        for inputs, _ in merged_test_cases:
            args = ["rgroup2smarts", "--single"]
            if check:
                args.append("--check")
            expected_output_lines = []
            for smiles in inputs:
                args.extend(["--cut-rgroup", smiles])
                output = simple_test_cases[smiles]
                expected_output_lines.append(output + "\n")
            stdout, stderr = run(args)
            expected_output = "".join(expected_output_lines)
            self.assertEqual(stdout, expected_output, repr(inputs))

    def test_single(self):
        self._test_single(False)

    def test_single_check(self):
        self._test_single(True)

    def _test_explain(self, check):
        args = ["rgroup2smarts", "--explain"]
        if check:
            args.append("--check")
        args.extend(
            [
                "--cut-rgroup",
                "c1ccccc1*",
                "--cut-rgroup",
                "O*",
            ]
        )
        stdout, stderr = run(args, True)
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1),$([OHv2])]\n")
        return stderr

    def test_explain(self):
        stderr = self._test_explain(check=False)
        self.assertEqual(
            stderr,
            """\
Using --cut-rgroup SMILES from the command-line
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
""",
        )

    def test_explain_check(self):
        stderr = self._test_explain(check=True)
        self.assertEqual(
            stderr,
            """\
Using --cut-rgroup SMILES from the command-line
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#1 passed the self-check
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
#2 passed the self-check
Checking that the SMARTS matches all of the input molecules
checked #0
checked #1
""",
        )


class TestSmilesFromFile(unittest.TestCase):
    def _test_merged(self, check):
        filename = support.create_test_filename(self, "groups.txt")

        for inputs, output in merged_test_cases:
            with open(filename, "w") as f:
                for smiles in inputs:
                    f.write(smiles + "\n")

            args = ["rgroup2smarts", filename]
            if check:
                args.append("--check")
            stdout, stderr = run(args)
            self.assertEqual(stdout, output + "\n", repr(inputs))

    def test_merged(self):
        self._test_merged(False)

    def test_merged_check(self):
        self._test_merged(True)

    def _test_single(self, check):
        filename = support.create_test_filename(self, "rgroups.txt")

        for inputs, _ in merged_test_cases:
            expected_output_lines = []
            with open(filename, "w") as f:
                for smiles in inputs:
                    f.write(smiles + "\n")
                    output = simple_test_cases[smiles]
                    expected_output_lines.append(output + "\n")
            expected_output = "".join(expected_output_lines)

            args = ["rgroup2smarts", "--single", filename]
            if check:
                args.append("--check")
            stdout, stderr = run(args)
            self.assertEqual(stdout, expected_output, repr(inputs))

    def test_single(self):
        self._test_single(False)

    def test_single_check(self):
        self._test_single(True)

    def _test_explain(self, check):
        filename = support.create_test_filename(self, "fra\tgs.smi")
        with open(filename, "w") as f:
            f.write("c1ccccc1*\nO*\n")
        args = ["rgroup2smarts", "--explain", filename]
        if check:
            args.append("--check")
        stdout, stderr = run(args, True)
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1),$([OHv2])]\n")

        return fix_stderr(self, stderr, filename)

    def test_explain(self):
        stderr = self._test_explain(check=False)
        self.assertEqual(
            stderr,
            """\
Reading R-group SMILES from frags.smi
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
""",
        )

    def test_explain_check(self):
        stderr = self._test_explain(check=True)
        self.assertEqual(
            stderr,
            """\
Reading R-group SMILES from frags.smi
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#1 passed the self-check
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
#2 passed the self-check
Checking that the SMARTS matches all of the input molecules
checked #0
checked #1
""",
        )

    def test_different_whitespace(self):
        filename = support.create_test_filename(self, "rgroups.txt")
        with open(filename, "w") as outfile:
            outfile.write("*Cl\tchlorine\n" "*Br bromine\n" "*F  and more\n")
        stdout, stderr = run(["rgroup2smarts", filename])
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, "*-!@[$([ClH0v1]),$([BrH0v1]),$([FH0v1])]\n")


# Basic burn-test that stdin works
class TestSmilesFromStdin(unittest.TestCase):
    def test_merged(self):
        stdout, stderr = run(["rgroup2smarts"], input="*c1ccccc1O\n*F\n")
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]\n")


######

_bad_smiles_inputs = [
    ("*Q", "Cannot parse SMILES ('*Q') at --cut-rgroup SMILES #1\n"),
    ("*C *Q", "Cannot parse SMILES ('*Q') at --cut-rgroup SMILES #2\n"),
    ("c1ccccc1", "Cannot convert SMILES ('c1ccccc1') at --cut-rgroup SMILES #1: no wildcard atom found\n"),
    ("*CN*", "Cannot convert SMILES ('*CN*') at --cut-rgroup SMILES #1: more than one wildcard atom\n"),
    (
        "*N[CH3:1]",
        "Cannot convert SMILES ('*N[CH3:1]') at --cut-rgroup SMILES #1: atom maps are not supported (atom 2 has atom map '1')\n",
    ),
    ("*", "Cannot convert SMILES ('*') at --cut-rgroup SMILES #1: wildcard atom not bonded to anything\n"),
    ("*N *=O", "Cannot convert SMILES ('*=O') at --cut-rgroup SMILES #2: wildcard atom not bonded via a single bond\n"),
    (
        "[*H]F",
        "Cannot convert SMILES ('[*H]F') at --cut-rgroup SMILES #1: wildcard atom must not have implicit hydrogens\n",
    ),
    ("[*-]F", "Cannot convert SMILES ('[*-]F') at --cut-rgroup SMILES #1: wildcard atom must be uncharged\n"),
    ("[*+2]F", "Cannot convert SMILES ('[*+2]F') at --cut-rgroup SMILES #1: wildcard atom must be uncharged\n"),
    ("Cl*F", "Cannot convert SMILES ('Cl*F') at --cut-rgroup SMILES #1: wildcard atom must only have one bond\n"),
    ("*Cl.F", "Cannot convert SMILES ('*Cl.F') at --cut-rgroup SMILES #1: more than one fragment found\n"),
]


class TestCommandlineFailures(unittest.TestCase):
    def test_bad_smiles(self):
        for smiles_list, errmsg in _bad_smiles_inputs:
            args = ["rgroup2smarts"]
            for smiles in smiles_list.split():
                args.extend(["--cut-rgroup", smiles])
            assert len(args) > 1, smiles_list
            stderr = run_fail(args)
            self.assertEqual(errmsg, stderr)


class TestFilenameFailures(unittest.TestCase):
    def test_bad_smiles(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        args = ["rgroup2smarts", filename]

        for smiles_list, errmsg in _bad_smiles_inputs:
            with open(filename, "w") as outfile:
                for smiles in smiles_list.split():
                    outfile.write(smiles + "\n")

            stderr = run_fail(args)
            stderr = fix_stderr(self, stderr, filename)
            errmsg = errmsg.replace("--cut-rgroup SMILES #", "frags.smi, line ")
            self.assertEqual(stderr, errmsg)

    def test_blank_line_not_allowed(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        with open(filename, "w") as f:
            f.write("*C\n\n*N\n")

        stderr = run_fail(["rgroup2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: no SMILES found at frags.smi, line 2\n")

    def test_blank_line_not_allowed(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        with open(filename, "w") as f:
            f.write("*C\n\n*N\n")

        stderr = run_fail(["rgroup2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: no SMILES found at frags.smi, line 2\n")

    def test_initial_whitespace_not_allowed(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        with open(filename, "w") as f:
            f.write("*C\n *N\n")

        stderr = run_fail(["rgroup2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: expected SMILES at start of line at frags.smi, line 2\n")

    def test_empty_file_not_allowed(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        with open(filename, "w") as f:
            f.close()

        stderr = run_fail(["rgroup2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot make a SMARTS: no SMILES strings found in frags.smi\n")

    def test_file_does_not_exist(self):
        filename = support.create_test_filename(self, "rgroups.dat")
        stderr = run_fail(["rgroup2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot open input file: [Errno 2] No such file or directory: frags.smi\n")


class TestOtherErrors(unittest.TestCase):
    def test_both_cut_rgroup_and_filename(self):
        stderr = run_fail(["rgroup2smarts", "rgroups.dat", "--cut-rgroup", "*F"])
        self.assertIn("Cannot specify both an R-group filename and a --cut-rgroup\n", stderr)


if __name__ == "__main__":
    unittest.main()
