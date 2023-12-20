from __future__ import print_function

import sys
import os
import tempfile
import shutil
from click.testing import CliRunner
from mmpdblib import cli


def expect_pass(args, input=None):
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, args, input=input)
    if result.exit_code:
        import shlex

        args_msg = " ".join(shlex.quote(word) for word in args)
        if result.exc_info:
            import traceback

            traceback.print_exception(*result.exc_info)
        raise AssertionError(f"SystemExit trying to run '{args_msg}': {result.exit_code}: {result.stderr}")
    return result


def expect_fail(args, input=None):
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, args, input=input)
    if not result.exit_code:
        raise AssertionError(f"Should have failed: {args!r}")
    return result


def get_filename(filename):
    return os.path.join(os.path.dirname(__file__), filename)


def create_testdir_and_filename(test_case, filename):
    dirname = tempfile.mkdtemp(prefix="mmpdb_test")
    test_case.addCleanup(shutil.rmtree, dirname)
    return dirname, os.path.join(dirname, filename)


def create_test_filename(test_case, filename):
    dirname = tempfile.mkdtemp(prefix="mmpdb_test")
    test_case.addCleanup(shutil.rmtree, dirname)
    return os.path.join(dirname, filename)
