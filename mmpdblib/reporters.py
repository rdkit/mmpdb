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

"""A 'reporter' is used to provide progress information and status reports"""

from __future__ import print_function, absolute_import

import sys
import time

from ._compat import basestring

def get_reporter(reporter):
    if reporter is None:
        return Quiet()
    if isinstance(reporter, basestring):
        if reporter == "quiet":
            return Quiet()
        if reporter == "verbose":
            return Verbose()
        raise ValueError("Unsupported reporter %r" % (reporter,))
    return reporter

class BaseReporter(object):
    def warning(self, msg):
        "Print a warning message"
        pass

    def report(self, msg):
        "Print a report message"
        pass

    def progress(self, it, text, n=None):
        "Return a context manager for giving status report about an iterator"
        return StatusContext(iter(it))

    def update(self, msg):
        "Update the status line. This will erase the previous status."
        pass

    def explain(self, msg, *args):
        if args:
            self.report(msg % args)
        else:
            self.report(msg)


class Quiet(BaseReporter):
    "This reporter does nothing"


# This lets me do things like:
class StatusContext(object):
    """Adapter to treat an iterator as a context manger"""
    def __init__(self, it):
        self._it = it

    def __enter__(self):
        return self

    def __iter__(self):
        return self._it

    def __exit__(self, *args):
        pass

class Verbose(BaseReporter):
    "This reporter sends report and status information to stderr."
    def __init__(self):
        self._erase = ""  # how to erase the last status message
        
    def warning(self, msg):
        "Clear any status message and report the warning"
        self.update("")
        sys.stderr.write("WARNING: %s\n" % (msg,))
        sys.stderr.flush()

    def report(self, msg):
        "Clear any status message and print the report line"
        if self._erase:
            self.update("")
        sys.stderr.write(msg + "\n")
        sys.stderr.flush()
        
        
    def progress(self, it, text, n=None):
        # Used in iterators

        def iterate():
            if n is None or n == 0:
                def get_text(i):
                    return text + " " + str(i)
            else:
                def get_text(i):
                    return text + " %d/%d (%.1f%%)" % (i, n, 100.0*i/n)

            i = 0
            self.update(get_text(i))
            t1 = time.time()
            try:
                for i, value in enumerate(it, 1):
                    yield value
                    t2 = time.time()
                    if t2 - t1 > 0.5:
                        self.update(get_text(i))
                        t1 = t2
            finally:
                self.update("")

        obj = StatusContext(iterate())
        # A bit of a hack so I can wrap location-based iterators
        if hasattr(it, "location"):
            obj.location = getattr(it, "location")
        return obj
            
    def update(self, msg):
        "Update the status line (erase the previous status message and display the new one)"
        sys.stderr.write(self._erase)
        sys.stderr.write(msg)
        sys.stderr.flush()
        self._erase = "\r" + " "*len(msg) + "\r"

# This is a bit of a hack that was developed at the very end of the project.
# It's used during the database load process, from an ".mmpa" file.

# It's an iter-like wrapper wihch shows progress across multiple
# stages.  It assumes that the total number of elements is known in
# the beginning, across all of the stages. For each stage, I want to
# show the overall progress as a percentage, and give some feedback
# about the progress of each stage.

class MultiStageReporter(object):
    def __init__(self, reporter, num_rows):
        self.reporter = reporter
        self.num_rows = num_rows
        self._it = None
        self._row_count = 0

    def set_iter(self, template, container):
        """A string template (must have the '%' terms in the right order) and the container to iterator over
        
        This must be called to start each stage.
        """
        self.template = template
        self._it = enumerate(container)  # enumerate() so I can track progress for the stage
        self._n = len(container)
        msg = self.template  % (100.0*self._row_count/self.num_rows, 0, self._n)
        self.reporter.update(msg)
        self._prev_time = time.time()
                 
    def __iter__(self):
        return self

    def __next__(self):
        try:
            i, value = next(self._it)
        except StopIteration:
            self.reporter.update("") # Reset at the end of the stage.
            raise
        self._row_count = row_count = self._row_count + 1  # Global count

        now = time.time()
        if now - self._prev_time > 0.5:
            # Show the progress. Template '%' terms must be: overall percentage,
            # element number in the stage, total number of elements in the stage.
            self.reporter.update(self.template % (100.0*row_count/self.num_rows, i, self._n))
            self._prev_time = now
        
        return value
        

    next = __next__
        
