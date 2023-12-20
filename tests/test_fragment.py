# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
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

import gzip
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import contextmanager


import mmpdblib
from mmpdblib import cli
from mmpdblib.config import DEFAULT_FRAGMENT_OPTIONS
from mmpdblib import fragment_db

from support import (
    get_filename,
    create_test_filename,
    expect_pass,
    expect_fail,
)

TAB_SMI = get_filename("tab.smi")
TWO_TABS_SMI = get_filename("two_tabs.smi")
SPACE_SMI = get_filename("space.smi")
SPACE_SMI_GZ = get_filename("space.smi.gz")
COMMA_SMI = get_filename("comma.smi")

# To make the cached.fragdb file, start with a SMILES file containing:
"""
Oc1ccccc1O record
Pc1ccccc1C entry
"""
# then process it with:
"""
python -m mmpdblib fragment cached.smi --max-heavies 22 --max-rotatable-bonds 40 \
     --rotatable-smarts '**' --cut-smarts '[R][\!R]' --num-cuts 2 \
     --salt-remover '<none>' -o cached.fragdb
"""
# (These options are used to test that the fragmentation options comes from the
# cache, and not from the default or command-line.)
#
# Then modify the Pc1ccccc1C so it appears to be Nc1ccccc1C
#
# % sqlite3 cached.fragdb
# sqlite> update record set input_smiles = "Nc1ccccc1C", normalized_smiles = "Cc1ccccc1N" where input_smiles = "Pc1ccccc1C";
# (This forces the cache lookup to be incorrect, which I use to test the cache is being used.)

CACHED_FRAGDB = get_filename("cached.fragdb")

from rdkit import Chem

wildcard_atom = Chem.CanonSmiles("*")
# Before 2018, RDKit returned "[*]"
assert wildcard_atom == "*", ("When did RDKit change?", wildcard_atom)


def working_dir(*filenames):
    dirname = tempfile.mkdtemp()
    try:
        for filename in filenames:
            basename = os.path.basename(filename)
            shutil.copy(filename, os.path.join(dirname, basename))
    except:
        shutil.rmtree(dirname)
        raise
    return WorkingDir(dirname)


class WorkingDir:
    def __init__(self, dirname):
        self.dirname = dirname

    def __enter__(self):
        return self

    def __exit__(self, *args):
        shutil.rmtree(self.dirname)

    def __truediv__(self, path):
        basename = os.path.basename(path)
        return os.path.join(self.dirname, basename)


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


# Handle automatic naming
_fragdb_name = {
    "tab.smi": "tab.fragdb",
    "comma.smi": "comma.fragdb",
    "space.smi": "space.fragdb",
    "two_tabs.smi": "two_tabs.fragdb",
    "space.smi.gz": "space.fragdb",
}


@contextmanager
def fragment(path, args, fragdb_name=None):
    if isinstance(args, str):
        args = args.split()
    assert isinstance(args, (tuple, list))
    with working_dir(path) as workdir:
        args = list(args) + [workdir / path]
        args = ("--quiet", "fragment") + fix_fragment_args(args)
        expect_pass(args, input=None)

        fragdb_name = _fragdb_name[os.path.basename(path)]
        with fragment_db.open_fragdb(workdir / fragdb_name) as db:
            yield db


@contextmanager
def fragment_stdin(text, args):
    if isinstance(args, str):
        args = args.split()
    assert isinstance(args, (tuple, list))
    with working_dir() as workdir:
        output_name = workdir / "testme.fragdb"
        args = list(args) + ["-o", output_name]
        args = ("--quiet", "fragment") + fix_fragment_args(args)
        expect_pass(args, input=text)

        with fragment_db.open_fragdb(output_name) as db:
            yield db


def fragment_fail(path, args):
    if isinstance(args, str):
        args = args.split()
    assert isinstance(args, (tuple, list))

    with working_dir(path) as workdir:
        args = ["--quiet", "fragment", "--salt-remover", "<none>"] + list(args) + [workdir / path]
        result = expect_fail(args)
        return result.stderr


def get_ids(records):
    return [record.id for record in records]


class TestSmilesParser(unittest.TestCase):

    # default, which is "whitespace"
    def test_space_as_default(self):
        with fragment(SPACE_SMI, []) as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_space_gz_as_default(self):
        with fragment(SPACE_SMI_GZ, []) as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_tab_as_default(self):
        with fragment(TAB_SMI, []) as db:
            self.assertEqual(get_ids(db), ["record", "entry"])

    def test_two_tabs_as_default(self):
        with fragment(TWO_TABS_SMI, []) as db:
            self.assertEqual(get_ids(db), ["record", "vinyl"])

    def test_comma_as_default(self):
        stderr = fragment_fail(COMMA_SMI, [])

        self.assertIn("must contain at least two whitespace-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)

    # "whitespace"
    def test_space_as_whitespace(self):
        with fragment(SPACE_SMI, "--delimiter whitespace") as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_space_gz_as_whitespace(self):
        with fragment(SPACE_SMI_GZ, "--delimiter whitespace") as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_tab_as_whitespace(self):
        with fragment(TAB_SMI, "--delimiter whitespace", "tab.fragdb") as db:
            self.assertEqual(get_ids(db), ["record", "entry"])

    def test_two_tabs_as_whitespace(self):
        with fragment(TWO_TABS_SMI, "--delimiter whitespace") as db:
            self.assertEqual(get_ids(db), ["record", "vinyl"])

    def test_comma_as_whitespace(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter whitespace")

        self.assertIn("must contain at least two whitespace-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)

    # "space"
    def test_space_as_space(self):
        with fragment(SPACE_SMI, "--delimiter space") as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_space_gz_as_space(self):
        with fragment(SPACE_SMI_GZ, "--delimiter space") as db:
            self.assertEqual(get_ids(db), ["record", "entry", "item"])

    def test_tab_as_space(self):
        with fragment(TAB_SMI, "--delimiter space") as db:
            self.assertEqual(get_ids(db), ["1", "2"])

    def test_two_tabs_as_space(self):
        stderr = fragment_fail(TWO_TABS_SMI, "--delimiter space")
        self.assertIn("must contain at least two space-delimited fields", stderr)

    def test_comma_as_space(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter space")
        self.assertIn("must contain at least two space-delimited fields", stderr)
        # The first line contains "Oc1ccccc1O,record 1", which is processed as whitespace
        # to give "Oc1ccccc1O,record" and "1". The first cannot be parsed as SMILES.
        # The second line has no space, which is why it has the error.
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)

    # "tab"
    def test_space_as_tab(self):
        stderr = fragment_fail(SPACE_SMI, ["--delimiter", "tab"])
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("space.smi", stderr)
        self.assertIn("line 1, record #1", stderr)

    def test_space_gz_as_tab(self):
        stderr = fragment_fail(SPACE_SMI_GZ, ["--delimiter", "tab"])
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("space.smi.gz", stderr)
        self.assertIn("line 1, record #1", stderr)

    def test_tab_as_tab(self):
        with fragment(TAB_SMI, "--delimiter tab") as db:
            self.assertEqual(get_ids(db), ["record 1", "entry 2"])

    def test_two_tabs_as_tab(self):
        with fragment(TWO_TABS_SMI, "--delimiter tab") as db:
            self.assertEqual(get_ids(db), ["record", "vinyl"])

    def test_comma_as_tab(self):
        stderr = fragment_fail(COMMA_SMI, "--delimiter tab")
        self.assertIn("must contain at least two tab-delimited fields", stderr)
        self.assertIn("Oc1ccccc1O,record 1", stderr)
        self.assertIn("line 1, record #1", stderr)

    # "to-eol"
    def test_space_as_to_eol(self):
        with fragment(SPACE_SMI, "--delimiter to-eol") as db:
            self.assertEqual(get_ids(db), ["record 1", "entry 2", "item 3"])

    def test_space_gz_as_to_eol(self):
        with fragment(SPACE_SMI_GZ, "--delimiter to-eol") as db:
            self.assertEqual(get_ids(db), ["record 1", "entry 2", "item 3"])

    def test_tab_as_to_eol(self):
        with fragment(TAB_SMI, "--delimiter to-eol") as db:
            self.assertEqual(get_ids(db), ["record 1", "entry 2"])

    def test_two_tabs_as_to_eol(self):
        with fragment(TWO_TABS_SMI, "--delimiter to-eol") as db:
            self.assertEqual(get_ids(db), ["record\t1", "vinyl\t2"])

    def test_comma_as_tab(self):
        stderr = fragment_fail(COMMA_SMI, ["--delimiter", "to-eol"])
        self.assertIn("must contain a whitespace to delimit the to-eol fields", stderr)
        self.assertIn("Nc1ccccc1C,entry,2", stderr)
        self.assertIn("line 2, record #2", stderr)

    # "comma"
    def test_space_as_comma(self):
        stderr = fragment_fail(SPACE_SMI, ["--delimiter", "comma"])
        self.assertIn("must contain at least two comma-delimited fields", stderr)

    def test_space_gz_as_comma(self):
        stderr = fragment_fail(SPACE_SMI_GZ, ["--delimiter", "comma"])
        self.assertIn("must contain at least two comma-delimited fields", stderr)

    def test_tab_as_comma(self):
        stderr = fragment_fail(TAB_SMI, ["--delimiter", "comma"])
        self.assertIn("must contain at least two comma-delimited fields", stderr)

    def test_two_tabs_as_comma(self):
        stderr = fragment_fail(TWO_TABS_SMI, ["--delimiter", "comma"])
        self.assertIn("must contain at least two comma-delimited fields", stderr)

    def test_comma_as_comma(self):
        with fragment(COMMA_SMI, "--delimiter comma") as db:
            self.assertEqual(get_ids(db), ["record 1", "entry", "item 3"])

    # --has-header
    def test_header_space(self):
        with fragment(SPACE_SMI, "--has-header") as db:
            self.assertEqual(get_ids(db), ["entry", "item"])

    def test_header_space_gz(self):
        with fragment(SPACE_SMI_GZ, "--has-header") as db:
            self.assertEqual(get_ids(db), ["entry", "item"])

    def test_header_tab(self):
        with fragment(TAB_SMI, "--delimiter tab --has-header") as db:
            self.assertEqual(get_ids(db), ["entry 2"])

    def test_header_comma(self):
        with fragment(COMMA_SMI, "--delimiter comma --has-header") as db:
            self.assertEqual(get_ids(db), ["entry", "item 3"])


class TestOptions(unittest.TestCase):
    def _check_error_record(self, rec, id, input_smiles, errmsg):
        self.assertEqual(rec.id, id)
        self.assertEqual(rec.input_smiles, input_smiles)
        self.assertEqual(rec.errmsg, errmsg)

    def _check_record(self, rec, id, input_smiles):
        self.assertEqual(rec.id, id)
        self.assertEqual(rec.input_smiles, input_smiles)
        self.assertEqual(rec.errmsg, None)

    ##  [--max-heavies N]
    def test_max_heavies(self):
        with fragment_stdin(
            "CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
            "--max-heavies 7",
        ) as db:
            self.assertEqual(db.options.max_heavies, 7)
            error_records = list(db.iter_error_records())
            self._check_error_record(error_records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many heavy atoms")
            records = list(db)
            self._check_record(records[0], "phenol", "c1ccccc1O")

        with fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n", ["--max-heavies", "6"]) as db:
            records = list(db.iter_error_records())
            self.assertEqual(len(records), 2)
            self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many heavy atoms")
            self._check_error_record(records[1], "phenol", "c1ccccc1O", "too many heavy atoms")

    def test_max_heavies_none(self):
        smiles = "C1" + "C" * 200 + "C1"
        smi = smiles + " R202\n"
        with fragment_stdin(smi, []) as db:
            records = list(db.iter_error_records())
            self._check_error_record(records[0], "R202", smiles, "too many heavy atoms")

        with fragment_stdin(smi, ["--max-heavies", "none"]) as db:
            records = list(db)
            self.assertIs(db.options.max_heavies, None)
            self._check_record(records[0], "R202", smiles)

    def test_max_heavies_errors(self):
        for term in ("-1", "A", "3.4"):
            stderr = fragment_fail(SPACE_SMI, ["--max-heavies", term])
            self.assertIn("Invalid value for '--max-heavies': must be a positive integer or 'none'", stderr)

    ##  [--max-rotatable-bonds N]
    def test_max_rotatable_bonds(self):
        with fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n", ["--max-rotatable-bonds", "30"]) as db:

            self.assertEqual(db.options.max_rotatable_bonds, 30)
            records = list(db)
            self.assertEqual(len(records), 2)
            self._check_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC")
            self._check_record(records[1], "phenol", "c1ccccc1O")

        with fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n", ["--max-rotatable-bonds", "6"]) as db:
            records = list(db.iter_error_records())
            self.assertEqual(len(records), 1)
            self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many rotatable bonds")

            records = list(db)
            self.assertEqual(len(records), 1)
            self._check_record(records[0], "phenol", "c1ccccc1O")

    def test_max_rotatable_bonds_none(self):
        smiles = "C" * 14
        smi = "C" * 14 + " C14\n"
        with fragment_stdin(smi, []) as db:
            records = list(db.iter_error_records())
            self._check_error_record(records[0], "C14", smiles, "too many rotatable bonds")

        with fragment_stdin(smi, ["--max-rotatable-bonds", "none"]) as db:
            self.assertIs(db.options.max_rotatable_bonds, None)
            records = list(db)
            self._check_record(records[0], "C14", smiles)

    def test_max_rotatable_bonds_errors(self):
        for term in ("-1", "A", "3.4"):
            stderr = fragment_fail(SPACE_SMI, ["--max-rotatable-bonds", term])
            self.assertIn("Invalid value for '--max-rotatable-bonds': must be a positive integer or 'none'", stderr)

    ##  [--rotatable-smarts SMARTS]
    def test_rotatable_smarts(self):
        with fragment_stdin("CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n", ["--max-rotatable-bonds", "3"]) as db:
            self.assertEqual(db.options.max_rotatable_bonds, 3)
            records = list(db)
            self.assertEqual(len(records), 1)
            self._check_record(records[0], "phenol", "c1ccccc1O")

            records = list(db.iter_error_records())
            self.assertEqual(len(records), 1)
            self._check_error_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC", "too many rotatable bonds")

        with fragment_stdin(
            "CCCCCCCCCCCCCCCCCCCC C20\nc1ccccc1O phenol\n",
            ["--max-rotatable-bonds", "3", "--rotatable-smarts", "[x1]-*"],  # only the end atoms
        ) as db:
            records = list(db)
            self.assertEqual(len(records), 2)
            self._check_record(records[0], "C20", "CCCCCCCCCCCCCCCCCCCC")
            self._check_record(records[1], "phenol", "c1ccccc1O")

    def test_rotatable_smarts_errors(self):
        BAD_SMARTS = "unable to parse SMARTS"
        BAD_ATOMS = "rotatable SMARTS must match exactly two atoms"
        BAD_BONDS = "rotatable SMARTS must connect both atoms"

        for term, errmsg in (("C", BAD_ATOMS), ("Q4", BAD_SMARTS), ("C.C", BAD_BONDS), ("***", BAD_ATOMS)):
            stderr = fragment_fail(SPACE_SMI, ["--rotatable-smarts", term])
            self.assertIn(errmsg, stderr)

    ## ##  [--salt-remover FILENAME]
    # TODO: write tests.
    ## def test_salt_remover(self):
    ##     pass

    ##  [--cut-smarts SMARTS]
    def test_cut_smarts(self):
        with fragment_stdin("CCCCCCC C7\n", ["--cut-smarts", "C-C", "--num-cuts", "1"]) as db:
            self.assertEqual(db.options.num_cuts, 1)
            records = list(db)
            self.assertEqual(len(records), 1)
            self._check_record(records[0], "C7", "CCCCCCC")
            self.assertEqual(len(records[0].fragmentations), 6)

        with fragment_stdin("CCCCCCC C7\n", ["--cut-smarts", "O-O"]) as db:
            records = list(db)
            self.assertEqual(len(records), 1)
            self._check_record(records[0], "C7", "CCCCCCC")
            self.assertEqual(len(records[0].fragmentations), 0)

    def test_cut_smarts_errors(self):
        BAD_SMARTS = "unable to parse SMARTS"
        BAD_ATOMS = "cut SMARTS must match exactly two atoms"
        BAD_BONDS = "cut SMARTS must connect both atoms"

        for term, errmsg in (("C", BAD_ATOMS), ("Q4", BAD_SMARTS), ("C.C", BAD_BONDS), ("***", BAD_ATOMS)):
            stderr = fragment_fail(SPACE_SMI, ["--cut-smarts", term])
            self.assertIn(errmsg, stderr)

    ##  [--num-cuts {1,2,3}]
    def test_num_cuts(self):
        with fragment_stdin("CCCCCCC C7\n", ["--cut-smarts", "CC", "--num-cuts", "1"]) as db:
            self.assertEqual(db.options.num_cuts, 1)
            records = list(db)
            self.assertEqual(len(records[0].fragmentations), 6)

        with fragment_stdin("CCCCCCC C7\n", ["--cut-smarts", "CC", "--num-cuts", "2"]) as db:
            self.assertEqual(db.options.num_cuts, 2)
            records = list(db)
            self.assertEqual(len(records[0].fragmentations), 15)

    def test_num_cuts_errors(self):
        for term in ("0", "-1", "A", "4", "1.0"):
            stderr = fragment_fail(SPACE_SMI, ["--num-cuts", term])
            self.assertIn("Invalid value for '--num-cuts': ", stderr)
            self.assertIn("is not one of '1', '2', '3'", stderr)

    ## [--cache SOURCE]
    def test_cache(self):
        self._test_cache(CACHED_FRAGDB)

    def _test_cache(self, cache_filename):
        with fragment(SPACE_SMI, ["--no-salt-remover", "--cache", cache_filename]) as db:
            options = db.options
            self.assertEqual(options.max_heavies, 22)
            self.assertEqual(options.max_rotatable_bonds, 40)
            self.assertEqual(options.rotatable_smarts, "**")
            self.assertEqual(options.cut_smarts, "[R][!R]")
            self.assertEqual(options.num_cuts, 2)
            self.assertEqual(options.method, "chiral")
            self.assertEqual(options.salt_remover, "<none>")

            records = list(db)
            self._check_record(records[1], "entry", "Nc1ccccc1C")
            # This is a fake record to check that the value came from the cache.
            self.assertIn("P", records[1].fragmentations[0].constant_smiles)

            # This doesn't come from the cache.
            self._check_record(records[2], "item", "Nc1cc(S)ccc1C")

    ## [--num-jobs N]
    # def test_num_jobs(self): # TODO: this is hard to automate
    def test_num_jobs_errors(self):
        for term in ("-1", "0", "", "N", "2.0"):
            stderr = fragment_fail(SPACE_SMI, ["--num-jobs", term])
            self.assertIn("must be a positive integer", stderr)

    ## [-i FORMAT]
    # The right test for this is to read a gzip file from stdin.
    # That's hard to automate.


def smifrag(*args):
    args = ("--quiet", "smifrag") + fix_fragment_args(args)
    result = expect_pass(args)
    lines = result.output.splitlines(False)
    del lines[:3]
    return [[term.strip() for term in line.split("|")] for line in lines]


def smifrag_fail(*args):
    args = ("--quiet", "smifrag") + fix_fragment_args(args)
    result = expect_fail(args)
    return result.stderr


class TestSmiFrag(unittest.TestCase):
    def test_phenol(self):
        result = smifrag("c1ccccc1O")
        self.assertEqual(
            result,
            [
                ["1", "N", "1", "1", "*O", "0", "6", "1", "*c1ccccc1", "c1ccccc1"],
                ["1", "N", "6", "1", "*c1ccccc1", "0", "1", "1", "*O", "O"],
            ],
        )

    def test_max_heavies(self):
        errmsg = smifrag_fail("c1ccccc1O", "--max-heavies", "5")
        self.assertIn("Cannot parse smiles: too many heavy atoms", errmsg)

    def test_rotatable_bonds(self):
        errmsg = smifrag_fail("c1ccccc1O", "--max-rotatable-bonds", "3", "--rotatable-smarts", "**")
        self.assertIn("Cannot parse smiles: too many rotatable bonds", errmsg)

    def test_cut_smarts_using_cut_AlkylChains(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "cut_AlkylChains")
        self.assertEqual(
            result,
            [
                ["1", "N", "3", "1", "*CCC", "0", "1", "1", "*C", "C"],
                ["1", "N", "1", "1", "*C", "0", "3", "1", "*CCC", "CCC"],
                ["2", "N", "1", "11", "*C*", "01", "3", "12", "*C.*CC", "-"],
                ["2", "N", "2", "11", "*CC*", "01", "2", "11", "*C.*C", "-"],
                ["1", "N", "2", "1", "*CC", "0", "2", "1", "*CC", "CC"],
            ],
        )

    def test_cut_smarts(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "[CH2][CH2]")
        self.assertEqual(result, [["1", "N", "2", "1", "*CC", "0", "2", "1", "*CC", "CC"]])

    def test_num_cuts(self):
        result = smifrag("CCCC")
        self.assertEqual(result, [])
        result = smifrag("CCCC", "--cut-smarts", "cut_AlkylChains", "--num-cuts", "1")
        self.assertEqual(
            result,
            [
                ["1", "N", "3", "1", "*CCC", "0", "1", "1", "*C", "C"],
                ["1", "N", "1", "1", "*C", "0", "3", "1", "*CCC", "CCC"],
                ["1", "N", "2", "1", "*CC", "0", "2", "1", "*CC", "CC"],
            ],
        )


class TestFragmentCutRGroups(unittest.TestCase):
    def test_one_cut_rgroup(self):
        with fragment(SPACE_SMI, ["--cut-rgroup", "Oc1ccccc1*"], "space.fragdb") as db:
            records = list(db)
            self.assertEqual(len(records), 3)

            record = records[0]
            self.assertEqual(record.normalized_smiles, "Oc1ccccc1O")
            self.assertEqual(len(record.fragmentations), 3)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*O", "*c1ccccc1O", "*c1ccccc1*")))

            record = records[1]
            self.assertEqual(record.normalized_smiles, "Cc1ccccc1N")
            self.assertEqual(len(record.fragmentations), 0)

            record = records[2]
            self.assertEqual(record.normalized_smiles, "Cc1ccc(S)cc1N")
            self.assertEqual(len(record.fragmentations), 0)

    def test_two_cut_rgroups(self):
        with fragment(SPACE_SMI, ["--cut-rgroup", "Oc1ccccc1*", "--cut-rgroup", "*c1ccccc1N"], "space.fragdb") as db:
            records = list(db)

            self.assertEqual(len(records), 3)

            record = records[0]
            self.assertEqual(record.normalized_smiles, "Oc1ccccc1O")
            self.assertEqual(len(record.fragmentations), 3)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*O", "*c1ccccc1O", "*c1ccccc1*")))

            record = records[1]
            self.assertEqual(record.normalized_smiles, "Cc1ccccc1N")
            self.assertEqual(len(record.fragmentations), 2)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*C", "*c1ccccc1N")))

            record = records[2]
            self.assertEqual(record.normalized_smiles, "Cc1ccc(S)cc1N")
            self.assertEqual(len(record.fragmentations), 0)

    def test_invalid_cut_rgroup(self):
        stderr = fragment_fail(SPACE_SMI, ["--cut-rgroup", "Oc1ccccc1*", "--cut-rgroup", "c1ccccc1N"])
        self.assertEqual(stderr, "Cannot convert SMILES ('c1ccccc1N') at --cut-rgroup #2: no wildcard atom found\n")

    ### R-groups from file

    def test_cut_rgroup_filename(self):
        filename = create_test_filename(self, "rgroups.txt")
        with open(filename, "w") as f:
            f.write("*O\n" "*N\n")

        with fragment(SPACE_SMI, ["--cut-rgroup-file", filename], "space.fragdb") as db:
            records = list(db)
            record = records[0]
            self.assertEqual(record.normalized_smiles, "Oc1ccccc1O")
            self.assertEqual(len(record.fragmentations), 3)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*O", "*c1ccccc1O", "*c1ccccc1*")))

            record = records[1]
            self.assertEqual(record.normalized_smiles, "Cc1ccccc1N")
            self.assertEqual(len(record.fragmentations), 2)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*N", "*c1ccccc1C")))

            record = records[2]
            self.assertEqual(record.normalized_smiles, "Cc1ccc(S)cc1N")
            self.assertEqual(len(record.fragmentations), 2)
            variable_smiles_set = set(frag.variable_smiles for frag in record.fragmentations)
            self.assertEqual(variable_smiles_set, set(("*N", "*c1cc(S)ccc1C")))

    def test_missing_rgroup_filename(self):
        filename = create_test_filename(self, "rgroups.txt")
        stderr = fragment_fail(SPACE_SMI, ["--cut-rgroup-file", filename])
        self.assertIn("Cannot use --cut-rgroup-file", stderr)
        self.assertIn(repr(filename), stderr)
        self.assertIn("No such file or directory", stderr)


class TestSmiFragCutRGroups(unittest.TestCase):
    def test_one_cut_rgroup(self):
        result = smifrag("Oc1ccccc1N", "--cut-rgroup", "Oc1ccccc1*")
        self.assertEqual(
            result,
            [
                ["1", "N", "1", "1", "*N", "0", "7", "1", "*c1ccccc1O", "Oc1ccccc1"],
                ["1", "N", "7", "1", "*c1ccccc1O", "0", "1", "1", "*N", "N"],
            ],
        )

    def test_two_cut_rgroups(self):
        result = smifrag("Cc1ccccc1N", "--cut-rgroup", "Oc1ccccc1*", "--cut-rgroup", "*c1ccccc1N")
        self.assertEqual(
            result,
            [
                ["1", "N", "7", "1", "*c1ccccc1N", "0", "1", "1", "*C", "C"],
                ["1", "N", "1", "1", "*C", "0", "7", "1", "*c1ccccc1N", "Nc1ccccc1"],
            ],
        )

    def test_invalid_cut_rgroup(self):
        errmsg = smifrag_fail("Cc1ccccc1N", "--cut-rgroup", "Oc1ccccc1*", "--cut-rgroup", "FCl")
        self.assertEqual(errmsg, "Cannot convert SMILES ('FCl') at --cut-rgroup #2: no wildcard atom found\n")

    ### R-groups from file

    def test_cut_rgroup_filename(self):
        filename = create_test_filename(self, "rgroups.txt")
        with open(filename, "w") as f:
            f.write("*O\n" "*N\n")

        result = smifrag("Cc1cc(O)ccc1N", "--cut-rgroup-file", filename)
        self.assertEqual(
            result,
            [
                ["1", "N", "1", "1", "*O", "0", "8", "1", "*c1ccc(N)c(C)c1", "Cc1ccccc1N"],
                ["1", "N", "8", "1", "*c1ccc(N)c(C)c1", "0", "1", "1", "*O", "O"],
                ["2", "N", "7", "12", "*c1ccc(*)c(C)c1", "10", "2", "12", "*N.*O", "-"],
                ["1", "N", "1", "1", "*N", "0", "8", "1", "*c1ccc(O)cc1C", "Cc1cccc(O)c1"],
                ["1", "N", "8", "1", "*c1ccc(O)cc1C", "0", "1", "1", "*N", "N"],
            ],
        )

    def test_missing_rgroup_filename(self):
        filename = create_test_filename(self, "rgroups.txt")
        stderr = smifrag_fail("Nc1ccccc1O", "--cut-rgroup-file", filename)
        self.assertIn("Cannot use --cut-rgroup-file", stderr)
        self.assertIn(repr(filename), stderr)
        self.assertIn("No such file or directory", stderr)


# TODO: Add a large number of tests for the different interesting ways to fragment.
#  - check for chiral
#  - check for directional bonds
#  - probably other things too
# Hmm, perhaps that should go into "test_fragment_algorithm.py"

if __name__ == "__main__":
    unittest.main()
