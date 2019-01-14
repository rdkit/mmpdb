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

from mmpdblib import commandline
from mmpdblib import dbutils

from support import get_filename, create_test_filename

TEST_DATA_FRAGMENTS = get_filename("test_data.fragments")
TEST_DATA_CSV = get_filename("test_data.csv")

def index(mmpdb_filename, *args):
    args = ("--quiet", "index", TEST_DATA_FRAGMENTS, "-o", mmpdb_filename) + tuple(args)
    try:
        commandline.main(args)
    except SystemExit as err:
        raise AssertionError("SystemExit trying to run %r: %s" % (args, err))
        

class TestIndexCommandline(unittest.TestCase):
    def _get_options(self, *args):
        mmpdb_filename = create_test_filename(self, "default.mmpdb")
        index(mmpdb_filename, *args)
        db = dbutils.open_database(mmpdb_filename)
        dataset = db.get_dataset()
        return dataset, json.loads(dataset.index_options_str)
    
    def test_default_filters(self):
        dataset, options = self._get_options()
        self.assertEqual(options,
                         {u'symmetric': False,
                          u'max_variable_heavies': 10})
        self.assertEqual(dataset.get_num_rules(), 47)
        self.assertEqual(dataset.get_num_pairs(), 342)
        self.assertEqual(dataset.get_property_names(), [])
        
    def test_set_filters(self):
        dataset, options = self._get_options(
            "--min-variable-heavies", "1",
            "--max-variable-heavies", "29",
            "--min-variable-ratio", "0.1",
            "--max-variable-ratio", "0.99",
            "--max-heavies-transf", "25",
            "--symmetric",
            "--max-frac-trans", "3")
        self.assertEqual(options, {
            u'symmetric': True,
            u'max_frac_trans': 3.0,
            u'max_heavies_transf': 25,
            u'max_variable_heavies': 29,
            u'max_variable_ratio': 0.99,
            u'min_variable_heavies': 1,
            u'min_variable_ratio': 0.1})
        self.assertEqual(dataset.get_num_rules(), 2*47) # because --symmetric
        self.assertEqual(dataset.get_num_pairs(), 2*342) 
        self.assertEqual(dataset.get_property_names(), [])

    def test_max_variable_heavies_none(self):
        dataset, options = self._get_options("--max-variable-heavies", "none")
        self.assertEqual(options, {
            u'symmetric': False,
            })
        
    def test_with_properties(self):
        dataset, options = self._get_options("--properties", TEST_DATA_CSV, "--title", "test data")
        self.assertEqual(dataset.title, "test data")
        self.assertEqual(dataset.get_property_names(), ["MW", "MP"])

    
        
if __name__ == "__main__":
    unittest.main()
