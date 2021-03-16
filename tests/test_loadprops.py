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


def loadprops(mmpdb_filename, *args):
    args = ("--quiet", "loadprops", "-p", TEST_DATA_CSV, mmpdb_filename) + tuple(args)
    try:
        commandline.main(args)
    except SystemExit as err:
        raise AssertionError("SystemExit trying to run %r: %s" % (args, err))


class TestLoadpropsCommandline(unittest.TestCase):
    def _get_options(self, *args):
        mmpdb_filename = create_test_filename(self, "default.mmpdb")
        index(mmpdb_filename, *args)
        db = dbutils.open_database(mmpdb_filename)
        dataset = db.get_dataset()
        return dataset, json.loads(dataset.index_options_str), mmpdb_filename

    def test_loadprops(self):
        dataset, options, mmpdb_filename = self._get_options()
        self.assertEqual(dataset.get_num_rule_environment_stats(), 0)
        loadprops(mmpdb_filename)
        self.assertEqual(dataset.get_num_rule_environment_stats(), 533)


if __name__ == "__main__":
    unittest.main()
