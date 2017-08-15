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

import sys

if sys.version_info.major == 2:
    # Python 2
    import __builtin__
    basestring = __builtin__.basestring
    
    def open_universal(filename):
        return open(filename, "rU")

    # Treat a gzip binary file as a text files.
    # This may cause problems mixing newline conventions.
    def io_wrapper(fileobj):
        return fileobj

    # gzip requires a binary file
    binary_stdin = sys.stdin
    binary_stdout = sys.stdout

    # Lazy map available through itertools
    import itertools
    imap = itertools.imap
    
else:
    # Python 3
    import io
    basestring = (str, bytes)

    def open_universal(filename):
        return open(filename, "r", newline=None)

    # Convert binary file into text file
    def io_wrapper(fileobj):
        return io.TextIOWrapper(fileobj, newline=None)

    # gzip requires a binary file
    binary_stdin = sys.stdin.buffer
    binary_stdout = sys.stdout.buffer
    
    # Lazy map is the default
    imap = map
    
