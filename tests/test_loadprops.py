import unittest
import json

from mmpdblib import dbutils

from support import (
    get_filename,
    create_test_filename,
    expect_pass,
)

TEST_DATA_FRAGDB = get_filename("test_data.fragdb")
TEST_DATA_CSV = get_filename("test_data.csv")


def index(mmpdb_filename, *args):
    args = ("--quiet", "index", TEST_DATA_FRAGDB, "-o", mmpdb_filename) + tuple(args)
    expect_pass(args)


def loadprops(mmpdb_filename, *args):
    args = ("--quiet", "loadprops", "-p", TEST_DATA_CSV, mmpdb_filename) + tuple(args)
    expect_pass(args)


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
