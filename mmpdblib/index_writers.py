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

from __future__ import absolute_import, print_function

# Some low-level/backend writers for MMPWriter.
# The options are:
#   - save as a simple line-oriented format
#   - save the data to a SQLite database using Python's sqlite3 library
#   - save the data to a SQLite database using the third-party apsw library.

import os
import sqlite3
import datetime
import json
import collections
import itertools
import tempfile
import atexit
import shutil

# To install apsw, do:
# pip install --user https://github.com/rogerbinns/apsw/releases/download/3.16.2-r1/apsw-3.16.2-r1.zip \
# --global-option=fetch --global-option=--version --global-option=3.16.2 --global-option=--all \
# --global-option=build --global-option=--enable-all-extensions

try:
    import apsw
except ImportError:
    apsw = None
    ## The speedup is only about 10%. Not enough to make a firm suggestion.
    ## import sys
    ## sys.stderr.write("You should install apsw.\n")
    import sqlite3

from . import schema
from .fragment_algorithm import get_num_heavies_from_smiles

nan = float("nan")

class TableIndexWriter(object):
    def __init__(self, outfile):
        self.outfile = outfile
        self._W = outfile.write

    def close(self):
        self.outfile.close()
        
    def rollback(self):
        self._W("ROLLBACK")
        self.close()
        
    def commit(self):
        self._W("COMMIT")
        self.close()
        
    def start(self, fragment_options, index_options):
        self._W("VERSION\tmmpa/3\n")
        self._W("FRAGMENT_OPTIONS\t%s\n" % (json.dumps(list(fragment_options.to_dict().items())),))
        self._W("INDEX_OPTIONS\t%s\n" % (json.dumps(list(index_options.to_dict().items())),))

    def add_property_name(self, property_name_idx, property_name):
        self._W("PROPNAME\t%d\t%s\n" % (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self._W("RULE_SMILES\t%d\t%s\n" % (smiles_idx, smiles))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self._W("RULE\t%d\t%d\t%d\n" % (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, environment_fingerprint):
        self._W("FINGERPRINT\t%d\t%s\n" % (fp_idx, environment_fingerprint))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self._W("RULEENV\t%d\t%d\t%d\t%d\n" % (rule_env_idx, rule_idx, env_fp_idx, radius))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self._W("COMPOUND\t%d\t%s\t%s\t%s\t%d\n" % (
            compound_idx, compound_id, input_smiles,
            normalized_smiles, num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self._W("CONSTANT_SMILES\t%d\t%s\n" % (smiles_idx, constant_smiles))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self._W("PAIR%d\t\t%d\t%d\t%d\t%d\n" % (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self._W("PROP\t%d\t%d\t%s\n" % (compound_idx, property_name_idx, value))

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        self._W("RULEENV_STATS\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                    ((rule_env_idx, property_name_idx) + tuple(values)))

    def end(self, reporter):
        pass



def open_table_index_writer(outfile):
    return TableIndexWriter(outfile)
    
        

class BaseSqliteIndexWriter(object):
    def __init__(self, db, conn, title):
        self.db = db
        self.conn = conn
        self.title = title

    def start(self, fragment_options, index_options):
        creation_date = datetime.datetime.now().isoformat(sep=" ")
        fragment_options_str = json.dumps(fragment_options.to_dict())
        index_options_str = json.dumps(index_options.to_dict())
    
        self.conn.execute("INSERT INTO dataset (id, mmpdb_version, title, creation_date, "
                          "    fragment_options, index_options, is_symmetric) "
                          "    VALUES (1, 1, ?, ?, ?, ?, ?)",
                        (self.title, creation_date, fragment_options_str,
                         index_options_str, index_options.symmetric))

    def add_property_name(self, property_name_idx, property_name):
        self.conn.execute("INSERT INTO property_name (id, name) VALUES (?, ?)",
                          (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self.conn.execute("INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES (?, ?, ?)",
                          (smiles_idx, smiles, get_num_heavies_from_smiles(smiles)))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self.conn.execute("INSERT INTO rule (id, from_smiles_id, to_smiles_id) "
                          "  VALUES (?, ?, ?)",
                          (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, environment_fingerprint):
        self.conn.execute("INSERT INTO environment_fingerprint (id, fingerprint) "
                          " VALUES (?, ?)",
                          (fp_idx, environment_fingerprint))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self.conn.execute("INSERT INTO rule_environment (id, rule_id, environment_fingerprint_id,  radius) "
                          "  VALUES (?, ?, ?, ?)",
                          (rule_env_idx, rule_idx, env_fp_idx, radius))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self.conn.execute("INSERT INTO compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
                          "   VALUES (?, ?, ?, ?, ?)",
                          (compound_idx, compound_id, input_smiles, normalized_smiles, num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self.conn.execute("INSERT INTO constant_smiles (id, smiles) VALUES (?, ?)",
                          (smiles_idx, constant_smiles))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self.conn.execute("INSERT INTO pair (id, rule_environment_id, compound1_id, compound2_id, constant_id) "
                          "  VALUES (?, ?, ?, ?, ?)",
                          (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self.conn.execute("INSERT INTO compound_property (compound_id, property_name_id, value) VALUES (?, ?, ?)",
                          (compound_idx, property_name_idx, value))

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = values
        assert rule_env_idx is not None
        assert property_name_idx is not None
        assert count is not None
        assert avg is not None
        assert min is not None
        assert q1 is not None
        assert median is not None
        assert q3 is not None
        assert max is not None
        # XXX check for too-large/infinite values?

        self.conn.execute("INSERT INTO rule_environment_statistics "
                          "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
                          "       skewness, min, q1, median, q3, max, paired_t, p_value) "
                          "  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                          (rule_env_idx, property_name_idx, count, avg, std,
                           kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value))
        
    def end(self, reporter):
        reporter.update("Building index ...")
        schema.create_index(self.conn)
        
        # Improve SQLite query planning
        reporter.update("Analyzing database ...")
        self.conn.execute("ANALYZE")
        
        reporter.update("Computing sizes ...")
        num_compounds = schema._get_one(self.conn.execute("SELECT count(*) from compound"))
        num_rules = schema._get_one(self.conn.execute("SELECT count(*) from rule"))
        num_pairs = schema._get_one(self.conn.execute("SELECT count(*) from pair"))
        num_envs = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment"))
        num_stats = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment_statistics"))
        self.conn.execute("UPDATE dataset set num_compounds=?, num_rules=?, num_pairs=?, "
                          "num_rule_environments=?, num_rule_environment_stats=? WHERE id = 1",
                          (num_compounds, num_rules, num_pairs, num_envs, num_stats))
        
        reporter.update("")

class SQLiteIndexWriter(BaseSqliteIndexWriter):
    def close(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def commit(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def rollback(self):
        self.conn.close()
        self.db.close()
    
class APSWIndexWriter(BaseSqliteIndexWriter):
    def start(self, fragment_options, index_options):
        self.conn.execute("BEGIN TRANSACTION")
        super(APSWIndexWriter, self).start(fragment_options, index_options)
    
    def close(self):
        self.conn.close()
        self.db.execute("COMMIT")
        self.db.close()

    def commit(self):
        self.conn.execute("COMMIT")
        self.conn.close()
        self.db.close()

    def rollback(self):
        #self.conn.execute("ROLLBACK")
        self.conn.close()
        self.db.close()
    
        
def open_sqlite_index_writer(filename, title):
    if filename != ":memory:":
        if os.path.exists(filename):
            os.unlink(filename)
    if apsw is None:
        db = sqlite3.connect(filename)
        klass = SQLiteIndexWriter
    else:
        db = apsw.Connection(filename)
        klass = APSWIndexWriter
    
    schema.create_schema_for_sqlite(db)
    conn = db.cursor()

    return klass(db, conn, title)

