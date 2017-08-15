from __future__ import print_function

import sys
import os
import tempfile
import shutil

if sys.version_info.major == 2:
    from io import BytesIO
    class StringIO(BytesIO):
        def write(self, s):
            if isinstance(s, unicode):
                s = s.encode("utf8")
            return super(StringIO, self).write(s)
else:
    from io import StringIO


class redirect_stdin(object):
    def __init__(self, text):
        self.text = text
    def __enter__(self):
        self._real_stdin = sys.stdin
        self.stream = sys.stdin = StringIO(self.text)
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdin = self._real_stdin
        self.stream.close()

class capture_stdout(object):
    def __enter__(self):
        self._real_stdout = sys.stdout
        self.stream = sys.stdout = StringIO()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self._real_stdout
        self.value = self.stream.getvalue()
        self.stream.close()

class capture_stderr(object):
    def __enter__(self):
        self._real_stderr = sys.stderr
        self.stream = sys.stderr = StringIO()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr = self._real_stderr
        self.value = self.stream.getvalue()
        self.stream.close()
        
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
