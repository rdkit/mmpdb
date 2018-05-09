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

import mmpdblib
from mmpdblib import commandline
from mmpdblib import fragment_io
from mmpdblib.config import DEFAULT_FRAGMENT_OPTIONS

from support import (
    redirect_stdin, capture_stdout, capture_stderr,
    get_filename, create_testdir_and_filename, StringIO)

TAB_SMI = get_filename("tab.smi")
TWO_TABS_SMI = get_filename("two_tabs.smi")
SPACE_SMI = get_filename("space.smi")
SPACE_SMI_GZ = get_filename("space.smi.gz")
COMMA_SMI = get_filename("comma.smi")
CACHED_FRAGMENTS = get_filename("cached.fragments")

from rdkit import Chem
wildcard_atom = Chem.CanonSmiles("*")
if wildcard_atom == "[*]":
    # The tests were made for the pre-2018 when RDKit returned "[*]"
    def FIX(s):
        return s
elif wildcard_atom == "*":
    # RDKit 2018 change the behavior to return "*"
    def FIX(s):
        return s.replace("[*]", "*")
else:
    raise AssertionError(wildcard_atom)

def fix_fragment_args(args):
    salt_remover = True
    new_args = []
    for arg in args:
        assert isinstance(arg, str), arg
        if arg == "--no-salt-remover":
            salt_remover = False
        else:
            new_args.append(arg)
            
    args = tuple(new_args)
    if salt_remover:
        args = ("--salt-remover", "<none>") + args
    return args

def fragment_to_stdout(*args):
    args = ("--quiet", "fragment") + fix_fragment_args(args)
    
    try:
        with capture_stdout() as stdout:
            commandline.main(args)
        
        return stdout.value
    except SystemExit as err:
        raise AssertionError("SystemExit trying to run %r: %s" % (args, err))

def fragment(*args):
    stdout = fragment_to_stdout(*args)
    try:
        reader = fragment_io.read_fragment_records(StringIO(stdout))
        return (reader, list(reader))
    except Exception as err:
        raise AssertionError("Error parsing results of %r: %s" % (args, err))

def fragment_stdin(text, *args):
    with redirect_stdin(text):
        return fragment(*args)
    
    
def fragment_fail(*args):
    try:
        with capture_stdout() as stdout, capture_stderr() as stderr:
            commandline.main(("--quiet", "fragment", "--salt-remover", "<none>") + args)
    except SystemExit:
        return stderr.value
    raise AssertionError("Should have failed: %r" % (args,))

    


class FragmentHeader(unittest.TestCase):
    def test_header(self):
        reader, records = fragment(SPACE_SMI)
        self.assertEqual(reader.version, "mmpdb-fragment/2")
        self.assertEqual(reader.software, "mmpdb-" + mmpdblib.__version__)

        options = reader.options
        self.assertEqual(options.max_heavies, DEFAULT_FRAGMENT_OPTIONS.max_heavies)
        self.assertEqual(options.max_rotatable_bonds, DEFAULT_FRAGMENT_OPTIONS.max_rotatable_bonds)
        self.assertEqual(options.rotatable_smarts, DEFAULT_FRAGMENT_OPTIONS.rotatable_smarts)
        self.assertEqual(options.cut_smarts, DEFAULT_FRAGMENT_OPTIONS.cut_smarts)
        self.assertEqual(options.num_cuts, DEFAULT_FRAGMENT_OPTIONS.num_cuts)
        self.assertEqual(options.method, DEFAULT_FRAGMENT_OPTIONS.method)
        self.assertEqual(options.salt_remover, "<none>")

def get_ids(records):
    return [record.id for record in records]
    
class TestSmilesParser(unittest.TestCase):
    # default, which is "whitespace"
    def test_space_as_default(self):
        reader, records = fragment(SPACE_SMI)
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_space_gz_as_default(self):
        reader, records = fragment(SPACE_SMI_GZ)
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_tab_as_default(self):
        reader, records = fragment(TAB_SMI)
        self.assertEqual(get_ids(records), ["record", "entry"])
        
    def test_two_tabs_as_default(self):
        reader, records = fragment(TWO_TABS_SMI)
        self.assertEqual(get_ids(records), ["record", "vinyl"])

    def test_comma_as_default(self):
        stderr = fragment_fail(COMMA_SMI)
        self.assertIn("must contain at least two whitespace-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)
        
    # "whitespace"
    def test_space_as_whitespace(self):
        reader, records = fragment(SPACE_SMI, "--delimiter", "whitespace")
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_space_gz_as_whitespace(self):
        reader, records = fragment(SPACE_SMI_GZ, "--delimiter", "whitespace")
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_tab_as_whitespace(self):
        reader, records = fragment(TAB_SMI, "--delimiter", "whitespace")
        self.assertEqual(get_ids(records), ["record", "entry"])
        
    def test_two_tabs_as_whitespace(self):
        reader, records = fragment(TWO_TABS_SMI, "--delimiter", "whitespace")
        self.assertEqual(get_ids(records), ["record", "vinyl"])

    def test_comma_as_whitespace(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter", "whitespace")
        self.assertIn("must contain at least two whitespace-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)
        

    # "space"
    def test_space_as_space(self):
        reader, records = fragment(SPACE_SMI, "--delimiter", "space")
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_space_gz_as_space(self):
        reader, records = fragment(SPACE_SMI_GZ, "--delimiter", "space")
        self.assertEqual(get_ids(records), ["record", "entry", "item"])
        
    def test_tab_as_space(self):
        reader, records = fragment(TAB_SMI, "--delimiter", "space")
        self.assertEqual(get_ids(records), ["1", "2"])
        
    def test_two_tabs_as_space(self):
        stderr = fragment_fail(TWO_TABS_SMI, "--delimiter", "space")
        self.assertIn("must contain at least two space-delimited fields", stderr)

    def test_comma_as_space(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter", "space")
        self.assertIn("must contain at least two space-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)
        
        
    # "tab"
    def test_space_as_tab(self):
        stderr = fragment_fail(SPACE_SMI, "--delimiter", "tab")
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("space.smi", stderr)
        self.assertIn("line 1, record #1", stderr)
        
    def test_space_gz_as_tab(self):
        stderr = fragment_fail(SPACE_SMI_GZ, "--delimiter", "tab")
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("space.smi.gz", stderr)
        self.assertIn("line 1, record #1", stderr)
        
    def test_tab_as_tab(self):
        reader, records = fragment(TAB_SMI, "--delimiter", "tab")
        self.assertEqual(get_ids(records), ["record 1", "entry 2"])
        
    def test_two_tabs_as_tab(self):
        reader, records = fragment(TWO_TABS_SMI, "--delimiter", "tab")
        self.assertEqual(get_ids(records), ["record", "vinyl"])

    def test_comma_as_tab(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter", "tab")
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("Oc1ccccc1O,record 1", stderr)
        self.assertIn("line 1, record #1", stderr)
        
    # "to-eol"
    def test_space_as_to_eol(self):
        reader, records = fragment(SPACE_SMI, "--delimiter", "to-eol")
        self.assertEqual(get_ids(records), ["record 1", "entry 2", "item 3"])
        
    def test_space_gz_as_to_eol(self):
        reader, records = fragment(SPACE_SMI_GZ, "--delimiter", "to-eol")
        self.assertEqual(get_ids(records), ["record 1", "entry 2", "item 3"])
        
    def test_tab_as_to_eol(self):
        reader, records = fragment(TAB_SMI, "--delimiter", "to-eol")
        self.assertEqual(get_ids(records), ["record 1", "entry 2"])
        
    def test_two_tabs_as_to_eol(self):
        reader, records = fragment(TWO_TABS_SMI, "--delimiter", "to-eol")
        self.assertEqual(get_ids(records), ["record\t1", "vinyl\t2"])
        
    def test_comma_as_tab(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter", "to-eol")
        self.assertIn("must contain a whitespace to delimit the to-eol fields", stderr)
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)


    # "comma"
    def test_space_as_comma(self):
        stderr = fragment_fail(SPACE_SMI, "--delimiter", "comma")
        self.assertIn("must contain at least two comma-delimited fields", stderr)
        
    def test_space_gz_as_comma(self):
        stderr = fragment_fail(SPACE_SMI_GZ, "--delimiter", "comma")
        self.assertIn("must contain at least two comma-delimited fields", stderr)
        
    def test_tab_as_comma(self):
        stderr = fragment_fail(TAB_SMI, "--delimiter", "comma")
        self.assertIn("must contain at least two comma-delimited fields", stderr)
        
    def test_two_tabs_as_comma(self):
        stderr = fragment_fail(TWO_TABS_SMI, "--delimiter", "comma")
        self.assertIn("must contain at least two comma-delimited fields", stderr)

    def test_comma_as_comma(self):
        reader, records = fragment(COMMA_SMI, "--delimiter", "comma")
        self.assertEqual(get_ids(records), ["record 1", "entry", "item 3"])

    # --has-header
    def test_header_space(self):
        reader, records = fragment(SPACE_SMI, "--has-header")
        self.assertEqual(get_ids(records), ["entry", "item"])
        
    def test_header_space_gz(self):
        reader, records = fragment(SPACE_SMI_GZ, "--has-header")
        self.assertEqual(get_ids(records), ["entry", "item"])
        
    def test_header_tab(self):
        reader, records = fragment(TAB_SMI, "--delimiter", "tab", "--has-header")
        self.assertEqual(get_ids(records), ["entry 2"])
        
    def test_header_comma(self):
        reader, records = fragment(COMMA_SMI, "--delimiter", "comma", "--has-header")
        self.assertEqual(get_ids(records), ["entry", "item 3"])

class TestOptions(unittest.TestCase):
    def _check_error_record(self, rec, id, input_smiles, errmsg):
        self.assertEqual(rec.id, id)
        self.assertEqual(rec.input_smiles, input_smiles)
        self.assertEqual(rec.errmsg, errmsg)
        self.assertEqual(rec.fragments, [])
        
    def _check_record(self, rec, id, input_smiles):
        self.assertEqual(rec.id, id)
        self.assertEqual(rec.input_smiles, input_smiles)
        self.assertEqual(rec.errmsg, None)

    ##  [--max-heavies N]
    def test_max_heavies(self):
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-heavies", "7")
        self.assertEqual(reader.options.max_heavies, 7)
        self.assertEqual(len(records), 2)
        self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many heavy atoms")
        self._check_record(records[1], "phenol", "c1ccccc1O")
        
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-heavies", "6")
        self.assertEqual(len(records), 2)
        self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many heavy atoms")
        self._check_error_record(records[1], "phenol", "c1ccccc1O", "too many heavy atoms")

    def test_max_heavies_none(self):
        smiles = "C1" + "C"*200 + "C1"
        reader, records = fragment_stdin(smiles + " R202\n")
        self._check_error_record(records[0], "R202", smiles, "too many heavy atoms")
        
        reader, records = fragment_stdin(smiles + " R202\n", "--max-heavies", "none")
        self.assertIs(reader.options.max_heavies, 'none')
        self._check_record(records[0], "R202", smiles)
        
    def test_max_heavies_errors(self):
        for term in ("-1", "A", "3.4"):
            stderr = fragment_fail(SPACE_SMI, "--max-heavies", term)
            self.assertIn("argument --max-heavies: must be a positive integer or 'none'", stderr)

    ##  [--max-rotatable-bonds N]
    def test_max_rotatable_bonds(self):
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-rotatable-bonds", "30")
        self.assertEqual(reader.options.max_rotatable_bonds, 30)
        self.assertEqual(len(records), 2)
        self._check_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC")
        self._check_record(records[1], "phenol", "c1ccccc1O")
        
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-rotatable-bonds", "6")
        self.assertEqual(len(records), 2)
        self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many rotatable bonds")
        self._check_record(records[1], "phenol", "c1ccccc1O")

    def test_max_rotatable_bonds_none(self):
        smiles = "C"*14
        reader, records = fragment_stdin(smiles + " C14\n")
        self._check_error_record(records[0], "C14", smiles, "too many rotatable bonds")
        
        reader, records = fragment_stdin(smiles + " C14\n", "--max-rotatable-bonds", "none")
        self.assertIs(reader.options.max_rotatable_bonds, 'none')
        self._check_record(records[0], "C14", smiles)
        
    def test_max_rotatable_bonds_errors(self):
        for term in ("-1", "A", "3.4"):
            stderr = fragment_fail(SPACE_SMI, "--max-rotatable-bonds", term)
            self.assertIn("argument --max-rotatable-bonds: must be a positive integer or 'none'", stderr)

    ##  [--rotatable-smarts SMARTS]
    def test_rotatable_smarts(self):
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-rotatable-bonds", "3")
        self.assertEqual(reader.options.max_rotatable_bonds, 3)
        self.assertEqual(len(records), 2)
        self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many rotatable bonds")
        self._check_record(records[1], "phenol", "c1ccccc1O")
        
        reader, records = fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
                                         "--max-rotatable-bonds", "3",
                                         "--rotatable-smarts", "[x1]-*") # only the end atoms
        self.assertEqual(len(records), 2)
        self._check_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC")
        self._check_record(records[1], "phenol", "c1ccccc1O")

    def test_rotatable_smarts_errors(self):
        BAD_SMARTS = "unable to parse SMARTS"
        BAD_ATOMS = "rotatable SMARTS must match exactly two atoms"
        BAD_BONDS = "rotatable SMARTS must connect both atoms"
        
        for term, errmsg in (("C", BAD_ATOMS),
                             ("Q4", BAD_SMARTS),
                             ("C.C", BAD_BONDS),
                             ("***", BAD_ATOMS)):
            stderr = fragment_fail(SPACE_SMI, "--rotatable-smarts", term)
            self.assertIn(errmsg, stderr)

    ## ##  [--salt-remover FILENAME]
    # TODO: write tests.
    ## def test_salt_remover(self):
    ##     pass

    ##  [--cut-smarts SMARTS]
    def test_cut_smarts(self):
        reader, records = fragment_stdin("CCCCCCC C7\n",
                                         "--cut-smarts", "C-C", "--num-cuts", "1")
        self.assertEqual(reader.options.num_cuts, 1)
        self.assertEqual(len(records), 1)
        self._check_record(records[0], "C7", "CCCCCCC")
        self.assertEqual(len(records[0].fragments), 6)
        
        reader, records = fragment_stdin("CCCCCCC C7\n",
                                         "--cut-smarts", "O-O")
        self.assertEqual(len(records), 1)
        self._check_record(records[0], "C7", "CCCCCCC")
        self.assertEqual(len(records[0].fragments), 0)

    def test_cut_smarts_errors(self):
        BAD_SMARTS = "unable to parse SMARTS"
        BAD_ATOMS = "cut SMARTS must match exactly two atoms"
        BAD_BONDS = "cut SMARTS must connect both atoms"
        
        for term, errmsg in (("C", BAD_ATOMS),
                             ("Q4", BAD_SMARTS),
                             ("C.C", BAD_BONDS),
                             ("***", BAD_ATOMS)):
            stderr = fragment_fail(SPACE_SMI, "--cut-smarts", term)
            self.assertIn(errmsg, stderr)
        

    ##  [--num-cuts {1,2,3}]
    def test_num_cuts(self):
        reader, records = fragment_stdin("CCCCCCC C7\n", "--cut-smarts", "CC", "--num-cuts", "1")
        self.assertEqual(reader.options.num_cuts, 1)
        self.assertEqual(len(records[0].fragments), 6)
        
        reader, records = fragment_stdin("CCCCCCC C7\n", "--cut-smarts", "CC", "--num-cuts", "2")
        self.assertEqual(reader.options.num_cuts, 2)
        self.assertEqual(len(records[0].fragments), 15)

    def test_num_cuts_errors(self):
        for term in ("0", "-1", "A", "4", "1.0"):
            stderr = fragment_fail(SPACE_SMI, "--num-cuts", term)
            self.assertIn("--num-cuts: must be '1', '2', or '3'", stderr)
        
    ## [--cache SOURCE]
    def test_cache(self):
        self._test_cache(CACHED_FRAGMENTS)

    def _test_cache(self, cache_filename):
        reader, records = fragment(SPACE_SMI, "--no-salt-remover", "--cache", cache_filename)
        options = reader.options
        self.assertEqual(options.max_heavies, 22)
        self.assertEqual(options.max_rotatable_bonds, 40)
        self.assertEqual(options.rotatable_smarts, "**")
        self.assertEqual(options.cut_smarts, "[R][!R]")
        self.assertEqual(options.num_cuts, 2)
        self.assertEqual(options.method, "chiral")
        self.assertEqual(options.salt_remover, "<none>")

        self._check_record(records[1], "entry", "Nc1ccccc1C")
        # This is a fake record to check that the value came from the cache.
        self.assertIn("P", records[1].fragments[0].constant_smiles)

        # This doesn't come from the cache.
        self._check_record(records[2], "item", "Nc1cc(S)ccc1C")

    def test_cache_gz(self):
        dirname, filename_gz = create_testdir_and_filename(self, "cached.fragments.gz")
        with gzip.open(filename_gz, "wb") as outfile:
            with open(CACHED_FRAGMENTS, "rb") as infile:
                outfile.write(infile.read())
        self._test_cache(filename_gz)
    
    ## [--num-jobs N]
    #def test_num_jobs(self): # TODO: this is hard to automate
    def test_num_jobs_errors(self):
        for term in ("-1", "0", "", "N", "2.0"):
            stderr = fragment_fail(SPACE_SMI, "--num-jobs", term)
            self.assertIn("must be a positive integer", stderr)

    ## [-i FORMAT]
    # The right test for this is to read a gzip file from stdin.
    # That's hard to automate.
    
    ## [--output FILENAME]
    def test_output_to_fragments(self):
        dirname, filename = create_testdir_and_filename(self, "output.fragments")
        fragment_to_stdout(SPACE_SMI, "--output", filename, "--delimiter", "to-eol")
        reader = fragment_io.read_fragment_records(filename)
        records = list(reader)
        self._check_record(records[0], "record 1", "Oc1ccccc1O")
        
    def test_output_to_fragments_gz(self):
        dirname, filename = create_testdir_and_filename(self, "output.fragments.gz")
        fragment_to_stdout(SPACE_SMI, "--output", filename, "--delimiter", "to-eol")
        with open(filename, "rb") as infile:
            text = infile.read(4)
            self.assertEqual(text, b"\x1f\x8b\x08\x08")
        reader = fragment_io.read_fragment_records(filename)
        records = list(reader)
        self._check_record(records[0], "record 1", "Oc1ccccc1O")
        
    # These are human-readable/diagnostic fingerprints.
    # TODO: They should be removed for the future.
    def test_output_to_fraginfo(self):
        dirname, filename = create_testdir_and_filename(self, "output.fraginfo")
        fragment_to_stdout(SPACE_SMI, "--output", filename, "--delimiter", "to-eol")
        with open(filename) as infile:
            content = infile.read()
        self.assertIn("#heavies", content)
        
    def test_output_to_fraginfo_gz(self):
        dirname, filename = create_testdir_and_filename(self, "output.fraginfo.gz")
        fragment_to_stdout(SPACE_SMI, "--output", filename, "--delimiter", "to-eol")
        with open(filename, "rb") as infile:
            text = infile.read(4)
            self.assertEqual(text, b"\x1f\x8b\x08\x08")
        with gzip.open(filename) as infile:
            content = infile.read()
        self.assertIn(b"#heavies", content)

    ## [--out FORMAT]
    # I need a test for --out fragments.gz sent to stdout. That's hard to automate.

    
    
def smifrag(*args):
    args = ("--quiet", "smifrag") + fix_fragment_args(args)
    try:
        with capture_stdout() as stdout:
            commandline.main(args)

        lines = stdout.value.splitlines(False)
        del lines[:3]
        return [[term.strip() for term in line.split("|")] for line in lines]
        #return stdout.value
    except SystemExit as err:
        raise AssertionError("SystemExit trying to run %r: %s" % (args, err))

def smifrag_fail(*args):
    args = ("--quiet", "smifrag") + fix_fragment_args(args)
    with capture_stderr() as stderr:
        try:
            commandline.main(args)
        except SystemExit as err:
            pass
        else:
            raise AssertionError("Should have failed: %r" % (args,))
    return stderr.value
    
        
class TestSmiFrag(unittest.TestCase):
    def test_phenol(self):
        result = smifrag("c1ccccc1O")
        self.assertEqual(result, [
            ['1', 'N', '1', '1', FIX('[*]O'), '0', '6', '1', FIX('[*]c1ccccc1'), 'c1ccccc1'],
            ['1', 'N', '6', '1', FIX('[*]c1ccccc1'), '0', '1', '1', FIX('[*]O'), 'O'],
            ])

    def test_max_heavies(self):
        errmsg = smifrag_fail("c1ccccc1O", "--max-heavies", "5")
        self.assertIn("Cannot parse --smiles: too many heavy atoms", errmsg)
    
    def test_rotatable_bonds(self):
        errmsg = smifrag_fail("c1ccccc1O", "--max-rotatable-bonds", "3",
                              "--rotatable-smarts", "**")
        self.assertIn("Cannot parse --smiles: too many rotatable bonds", errmsg)
        
    def test_cut_smarts_using_cut_AlkylChains(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "cut_AlkylChains")
        self.assertEqual(result, [
            ['1', 'N', '3', '1', FIX('[*]CCC'), '0', '1', '1', FIX('[*]C'), 'C'],
            ['1', 'N', '1', '1', FIX('[*]C'), '0', '3', '1', FIX('[*]CCC'), 'CCC'],
            ['2', 'N', '1', '11', FIX('[*]C[*]'), '01', '3', '12', FIX('[*]C.[*]CC'), '-'],
            ['2', 'N', '2', '11', FIX('[*]CC[*]'), '01', '2', '11', FIX('[*]C.[*]C'), '-'],
            ['1', 'N', '2', '1', FIX('[*]CC'), '0', '2', '1', FIX('[*]CC'), 'CC']
            ])

    def test_cut_smarts(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "[CH2][CH2]")
        self.assertEqual(result, [
            ['1', 'N', '2', '1', FIX('[*]CC'), '0', '2', '1', FIX('[*]CC'), 'CC']
            ])

    def test_num_cuts(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "cut_AlkylChains",
                         "--num-cuts", "1")
        self.assertEqual(result, [
            ['1', 'N', '3', '1', FIX('[*]CCC'), '0', '1', '1', FIX('[*]C'), 'C'],
            ['1', 'N', '1', '1', FIX('[*]C'), '0', '3', '1', FIX('[*]CCC'), 'CCC'],
            ['1', 'N', '2', '1', FIX('[*]CC'), '0', '2', '1', FIX('[*]CC'), 'CC']
            ])
        
# TODO: Add a large number of tests for the different interesting ways to fragment.
#  - check for chiral
#  - check for directional bonds
#  - probably other things too
# Hmm, perhaps that should go into "test_fragment_algorithm.py"

if __name__ == "__main__":
    unittest.main()
