import sys
import unittest

from mmpdblib import commandline
from mmpdblib import fragment2smarts

import support

def run(cmd, stderr_ok=False):
    if isinstance(cmd, str):
        cmd = cmd.split()

    try:
        with support.capture_stdout() as stdout:
            with support.capture_stderr() as stderr:
                    commandline.main(cmd)
    except SystemExit as sys_exit:
        raise AssertionError(("unexpected exit", sys_exit, stderr.value))

    if not stderr_ok and stderr.value:
        print(stderr.value, file=sys.stderr)
        raise AssertionError(("unexpected stderr", cmd, stderr.value))
    return stdout.value, stderr.value

def run_failure(cmd):
    if isinstance(cmd, str):
        cmd = cmd.split()

    try:
        with support.capture_stdout() as stdout:
            with support.capture_stderr() as stderr:
                    commandline.main(cmd)
    except SystemExit as sys_exit:
        pass
    else:
        raise AssertionError("Failed to fail", cmd)

    return stdout.value, stderr.value
    

class TestSmilesOnCommandline(unittest.TestCase):
    def test_one_smiles(self):
        stdout, stderr = run("frag2smarts --cut-fragment *c1ccccc1O")
        self.assertEqual(
            stdout,
            '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2])]\n'
            )
    def test_two_smiles(self):
        stdout, stderr = run("frag2smarts --cut-fragment *c1ccccc1O --cut-fragment *F")
        self.assertEqual(
            stdout,
            '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]\n'
            )
    def test_four_smiles(self):
        stdout, stderr = run("frag2smarts --cut-fragment *c1ccccc1O --cut-fragment *F "
                                 "--cut-fragment *Cl --cut-fragment *[OH]")
        self.assertEqual(
            stdout,
            '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1]),$([ClH0v1]),$([OHv2])]\n'
            )
    def test_smiles_with_isotope_and_charge(self):
        stdout, stderr = run("frag2smarts --cut-fragment *c1ccccc1[16O] --cut-fragment *-[C+](=O)[O-]")
        self.assertEqual(
            stdout,
            '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]),$([C+H0v4](=[OH0v2])-[O-H0v1])]\n'
            )

def fix_stderr(test_case, stderr, filename):
    first_line, mid, rest = stderr.partition("\n")
    test_case.assertIn(repr(filename), first_line)
    first_line = first_line.replace(repr(filename), "frags.smi")

    return first_line + mid + rest
        

merged_test_cases = (
    (["*c1ccccc1O"], '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2])]'),
    (["*c1ccccc1O", "*F"], '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]'),
    (["*c1ccccc1O", "*F", "*Cl", "*[OH]"],
         '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1]),$([ClH0v1]),$([OHv2])]'),
    (["*c1ccccc1[16O]", "*-[C+](=O)[O-]"],
         '*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]),$([C+H0v4](=[OH0v2])-[O-H0v1])]'),
         )
    
simple_test_cases = {
    "*c1ccccc1O": '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]',
    "*F": '*-!@[FH0v1]',
    "*Cl": '*-!@[ClH0v1]',
    "*[OH]": '*-!@[OHv2]',
    "*c1ccccc1[16O]": '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]',
    "*-[C+](=O)[O-]": '*-!@[C+H0v4](=[OH0v2])-[O-H0v1]',
         }
    
        
class TestSmilesOnCommandline(unittest.TestCase):
    def _test_merged(self, check):
        for inputs, output in merged_test_cases:
            args = ["frag2smarts"]
            if check:
                args.append("--check")
            for smiles in inputs:
                args.extend(["--cut-fragment", smiles])
            stdout, stderr = run(args)
            self.assertEqual(stdout, output+"\n", repr(inputs))

    def test_merged(self):
        self._test_merged(False)
        
    def test_merged_check(self):
        self._test_merged(True)
    
    def _test_single(self, check):
        for inputs, _ in merged_test_cases:
            args = ["frag2smarts", "--single"]
            if check:
                args.append("--check")
            expected_output_lines = []
            for smiles in inputs:
                args.extend(["--cut-fragment", smiles])
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
        args = ["frag2smarts", "--explain"]
        if check:
            args.append("--check")
        args.extend([
                    "--cut-fragment", "c1ccccc1*",
                    "--cut-fragment", "O*",
            ])
        stdout, stderr = run(args, True)
        self.assertEqual(stdout,
                         "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1),$([OHv2])]\n")
        return stderr
    
    def test_explain(self):
        stderr = self._test_explain(check=False)
        self.assertEqual(stderr, """\
Using --cut-fragment SMILES from the command-line
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
""")

    def test_explain_check(self):
        stderr = self._test_explain(check=True)
        self.assertEqual(stderr, """\
Using --cut-fragment SMILES from the command-line
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#1 passed the self-check
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
#2 passed the self-check
Checking that the SMARTS matches all of the input molecules
checked #0
checked #1
""")


class TestSmilesFromFile(unittest.TestCase):
    def _test_merged(self, check):
        filename = support.create_test_filename(self, "fragments.txt")
        
        for inputs, output in merged_test_cases:
            with open(filename, "w") as f:
                for smiles in inputs:
                    f.write(smiles + "\n")
            
            args = ["frag2smarts", filename]
            if check:
                args.append("--check")
            stdout, stderr = run(args)
            self.assertEqual(stdout, output+"\n", repr(inputs))

    def test_merged(self):
        self._test_merged(False)
        
    def test_merged_check(self):
        self._test_merged(True)
        
    def _test_single(self, check):
        filename = support.create_test_filename(self, "fragments.txt")
        
        for inputs, _ in merged_test_cases:
            expected_output_lines = []
            with open(filename, "w") as f:
                for smiles in inputs:
                    f.write(smiles + "\n")
                    output = simple_test_cases[smiles]
                    expected_output_lines.append(output + "\n")
            expected_output = "".join(expected_output_lines)
                    
            args = ["frag2smarts", "--single", filename]
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
        args = ["frag2smarts", "--explain", filename]
        if check:
            args.append("--check")
        stdout, stderr = run(args, True)
        self.assertEqual(stdout,
                         "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1),$([OHv2])]\n")

        return fix_stderr(self, stderr, filename)
    
    def test_explain(self):
        stderr = self._test_explain(check=False)
        self.assertEqual(stderr, """\
Reading SMILES from frags.smi
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
""")

    def test_explain_check(self):
        stderr = self._test_explain(check=True)
        self.assertEqual(stderr, """\
Reading SMILES from frags.smi
#1: converted SMILES 'c1ccccc1*' to SMARTS '*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cHv4]:1'
#1 passed the self-check
#2: converted SMILES 'O*' to SMARTS '*-!@[OHv2]'
#2 passed the self-check
Checking that the SMARTS matches all of the input molecules
checked #0
checked #1
""")

# Basic burn-test that stdin works
class TestSmilesFromStdin(unittest.TestCase):
    def test_merged(self):
        with support.redirect_stdin("*c1ccccc1O\n*F\n"):
            stdout, stderr = run(["frag2smarts"])
        self.assertEqual(stdout, "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]\n")

######

_bad_smiles_inputs = [
    ("*Q", "Cannot parse SMILES ('*Q') at --cut-fragment SMILES #1\n"),
    ("*C *Q", "Cannot parse SMILES ('*Q') at --cut-fragment SMILES #2\n"),
    ("c1ccccc1", "Cannot convert SMILES ('c1ccccc1') at --cut-fragment SMILES #1: no wildcard atom found\n"),
    ("*CN*", "Cannot convert SMILES ('*CN*') at --cut-fragment SMILES #1: more than one wildcard atom\n"),
    ("*N[CH3:1]", "Cannot convert SMILES ('*N[CH3:1]') at --cut-fragment SMILES #1: atom maps are not supported (atom 2 has atom map '1')\n"),
    ("*.Cl", "Cannot convert SMILES ('*.Cl') at --cut-fragment SMILES #1: wildcard atom not bonded to anything\n"),
    ("*N *=O", "Cannot convert SMILES ('*=O') at --cut-fragment SMILES #2: wildcard atom not bonded via a single bond\n"),
    ("[*H]F", "Cannot convert SMILES ('[*H]F') at --cut-fragment SMILES #1: wildcard atom must not have implicit hydrogens\n"),
    ("[*-]F", "Cannot convert SMILES ('[*-]F') at --cut-fragment SMILES #1: wildcard atom must be uncharged\n"),
    ("[*+2]F", "Cannot convert SMILES ('[*+2]F') at --cut-fragment SMILES #1: wildcard atom must be uncharged\n"),
    ("Cl*F", "Cannot convert SMILES ('Cl*F') at --cut-fragment SMILES #1: wildcard atom must only have one bond\n"),
    ]

class TestCommandlineFailures(unittest.TestCase):
    def test_bad_smiles(self):
        for smiles_list, errmsg in _bad_smiles_inputs:
            args = ["frag2smarts"]
            for smiles in smiles_list.split():
                args.extend(["--cut-fragment", smiles])
            assert len(args) > 1, smiles_list
            stdout, stderr = run_failure(args)
            self.assertEqual(errmsg, stderr)

class TestFilenameFailures(unittest.TestCase):
    def test_bad_smiles(self):
        filename = support.create_test_filename(self, "fragments.dat")
        args = ["frag2smarts", filename]
        
        for smiles_list, errmsg in _bad_smiles_inputs:
            with open(filename, "w") as outfile:
                for smiles in smiles_list.split():
                    outfile.write(smiles + "\n")

            stdout, stderr = run_failure(args)
            stderr = fix_stderr(self, stderr, filename)
            errmsg = errmsg.replace("--cut-fragment SMILES #", "frags.smi, line ")
            self.assertEqual(stderr, errmsg)

    def test_blank_line_not_allowed(self):
        filename = support.create_test_filename(self, "fragments.dat")
        with open(filename, "w") as f:
            f.write("*C\n\n*N\n")

        stdout, stderr = run_failure(["frag2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: no SMILES found at frags.smi, line 2\n")

    def test_blank_line_not_allowed(self):
        filename = support.create_test_filename(self, "fragments.dat")
        with open(filename, "w") as f:
            f.write("*C\n\n*N\n")

        stdout, stderr = run_failure(["frag2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: no SMILES found at frags.smi, line 2\n")

    def test_initial_whitespace_not_allowed(self):
        filename = support.create_test_filename(self, "fragments.dat")
        with open(filename, "w") as f:
            f.write("*C\n *N\n")

        stdout, stderr = run_failure(["frag2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot parse input file: expected SMILES at start of line at frags.smi, line 2\n")

    def test_empty_file_not_allowed(self):
        filename = support.create_test_filename(self, "fragments.dat")
        with open(filename, "w") as f:
            f.close()

        stdout, stderr = run_failure(["frag2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot make a SMARTS: no SMILES strings found in frags.smi\n")

    def test_file_does_not_exist(self):
        filename = support.create_test_filename(self, "fragments.dat")
        stdout, stderr = run_failure(["frag2smarts", filename])
        stderr = fix_stderr(self, stderr, filename)
        self.assertEqual(stderr, "Cannot open input file: [Errno 2] No such file or directory: frags.smi\n")

class TestOtherErrors(unittest.TestCase):
    def test_both_cut_fragment_and_filename(self):
        stdout, stderr = run_failure(["frag2smarts", "fragments.dat", "--cut-fragment", "*F"])
        self.assertIn("Cannot specify both a fragment filename and a --cut-fragment\n", stderr)

if __name__ == "__main__":
    unittest.main()
    
