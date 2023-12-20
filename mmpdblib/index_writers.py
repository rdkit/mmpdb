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
import re


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
from . import fileio
from .fragment_algorithm import get_num_heavies_from_smiles

##nan = float("nan")

from . import index_algorithm

class CSVPairWriter(index_algorithm.BaseWriter):
    def start(self):
        self.num_pairs = 0
    
    def write_matched_molecule_pairs(self, pairs):
        backend = self.backend
        fragment_index = self.fragment_index

        n = 0
        for pair in pairs:
            rec1 = fragment_index.get_input_record(pair.id1)
            rec2 = fragment_index.get_input_record(pair.id2)
            backend.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                rec1.input_smiles, rec2.input_smiles, pair.id1, pair.id2, pair.smirks, pair.constant_smiles))
            n += 1
        self.num_pairs += n

    def __exit__(self, type, value, traceback):
        self.backend.close()


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
        self._W("VERSION\tmmpa/4\n")
        self._W("FRAGMENT_OPTIONS\t%s\n" % (json.dumps(list(fragment_options.to_dict().items())),))
        self._W("INDEX_OPTIONS\t%s\n" % (json.dumps(list(index_options.to_dict().items())),))

    def add_property_name(self, property_name_idx, property_name):
        self._W("PROPNAME\t%d\t%s\n" % (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self._W("RULE_SMILES\t%d\t%s\n" % (smiles_idx, smiles))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self._W("RULE\t%d\t%d\t%d\n" % (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, smarts, pseudosmiles, parent_smarts):
        if parent_smarts is None:
            parent_smarts = ""
        self._W("FINGERPRINT\t%d\t%s\t%s\t%s\n" % (fp_idx, smarts,  pseudosmiles, parent_smarts))

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

ADD_DATASET_SQL = (
    "INSERT INTO dataset (id, mmpdb_version, title, creation_date, "
    "    fragment_options, index_options, is_symmetric) "
    "    VALUES (1, 4, ?, ?, ?, ?, ?)"
    )
UPDATE_DATASET_SQL = (
    "UPDATE dataset SET \n"
    "    num_compounds=(SELECT COUNT(*) FROM compound), \n"
    "    num_rules=(SELECT COUNT(*) FROM rule), \n"
    "    num_pairs=(SELECT COUNT(*) FROM pair), \n"
    "    num_rule_environments=(SELECT COUNT(*) FROM rule_environment), \n"
    "    num_rule_environment_stats=(SELECT COUNT(*) FROM rule_environment_statistics) \n"
    "     WHERE id = 1"
    )

UPDATE_RULE_ENV_NUM_PAIRS = (
    "UPDATE rule_environment \n"
    "    SET num_pairs = ( \n"
    "       SELECT COUNT(*) FROM pair WHERE pair.rule_environment_id = rule_environment.id) \n"
    )
    
ADD_PROPERTY_NAME_SQL = (
    "INSERT INTO property_name (id, name) VALUES (?, ?)"
    )
ADD_RULE_SMILES_SQL = (
    "INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES (?, ?, ?)"
    )
ADD_RULE_SQL = (
    "INSERT INTO rule (id, from_smiles_id, to_smiles_id) VALUES (?, ?, ?)"
    )
ADD_ENVIRONMENT_FINGERPRINT_SQL = (
    "INSERT INTO environment_fingerprint (id, smarts, pseudosmiles, parent_smarts) VALUES (?, ?, ?, ?)"
    )
# num_pairs will always be set to zero but don't set it in the INSERT statement
# so it can be used by execute_values() in the pgsql adapter.
ADD_RULE_ENVIRONMENT_SQL = (
    "INSERT INTO "
    "  rule_environment (id, rule_id, environment_fingerprint_id,  radius, num_pairs) "
    "  VALUES (?, ?, ?, ?, ?)"
    )
ADD_COMPOUND_SQL = (
    "INSERT INTO "
    "  compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies)"
    "  VALUES (?, ?, ?, ?, ?)"
    )
ADD_CONSTANT_SMILES_SQL = (
    "INSERT INTO constant_smiles (id, smiles) VALUES (?, ?)"
    )
ADD_RULE_ENVIRONMENT_PAIR_SQL = (
    "INSERT INTO "
    "  pair (id, rule_environment_id, compound1_id, compound2_id, constant_id)"
    "  VALUES (?, ?, ?, ?, ?)"
    )
ADD_COMPOUND_PROPERTY_SQL = (
    "INSERT INTO compound_property (compound_id, property_name_id, value)"
    "  VALUES (?, ?, ?)"
    )
ADD_RULE_ENVIRONMENT_STATISTICS_SQL = (
    "INSERT INTO rule_environment_statistics "
    "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
    "    skewness, min, q1, median, q3, max, paired_t, p_value) "
    "  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    )

class StartMixin(object):
    def start(self, fragment_options, index_options):
        creation_date = datetime.datetime.now().isoformat(sep=" ")
        fragment_options_str = json.dumps(fragment_options.to_dict())
        index_options_str = json.dumps(index_options.to_dict())
        self.add_dataset(self.title, creation_date, fragment_options_str,
                        index_options_str, int(index_options.symmetric))
        
class BaseRDBMSIndexWriter(StartMixin):
    """Base class to save the index directly to a relational database management system"""
    def __init__(self, uri, db, conn, title):
        self.uri = uri
        self.db = db
        self.conn = conn
        self.title = title

        self.prepare_sql_statements()

    def prepare_sql_statements(self):
        # In theory this could actually execute prepare statments.
        # In practice, I use it to convert "?" to "%s" as appropriate
        self.ADD_DATASET_SQL = self.prepare_single(ADD_DATASET_SQL)
        
        self.ADD_PROPERTY_NAME_SQL = self.prepare_many(ADD_PROPERTY_NAME_SQL)
        self.ADD_RULE_SMILES_SQL = self.prepare_many(ADD_RULE_SMILES_SQL)
        self.ADD_RULE_SQL = self.prepare_many(ADD_RULE_SQL)
        self.ADD_ENVIRONMENT_FINGERPRINT_SQL = self.prepare_many(ADD_ENVIRONMENT_FINGERPRINT_SQL)
        self.ADD_RULE_ENVIRONMENT_SQL = self.prepare_many(ADD_RULE_ENVIRONMENT_SQL)
        self.ADD_COMPOUND_SQL = self.prepare_many(ADD_COMPOUND_SQL)
        self.ADD_CONSTANT_SMILES_SQL = self.prepare_many(ADD_CONSTANT_SMILES_SQL)
        self.ADD_RULE_ENVIRONMENT_PAIR_SQL = self.prepare_many(ADD_RULE_ENVIRONMENT_PAIR_SQL)
        self.ADD_COMPOUND_PROPERTY_SQL = self.prepare_many(ADD_COMPOUND_PROPERTY_SQL)
        self.ADD_RULE_ENVIRONMENT_STATISTICS_SQL = self.prepare_many(ADD_RULE_ENVIRONMENT_STATISTICS_SQL)

    def prepare_single(self, sql):
        # These will only be used in "execute()"
        return self.prepare(sql)

    def prepare_many(self, sql):
        # These will only be used in "execute_many()"
        return self.prepare(sql)

    def prepare(self, sql):
        return sql
        

    def create_schema(self):
        raise NotImplementedError

    ## def _execute(self, sql, args=()):
    ##     self.conn.execute(sql, args)
    ##     return self.conn

    ## def start(self, fragment_options, index_options):
    ##     creation_date = datetime.datetime.now().isoformat(sep=" ")
    ##     fragment_options_str = json.dumps(fragment_options.to_dict())
    ##     index_options_str = json.dumps(index_options.to_dict())
    ##     self.add_dataset(self.title, creation_date, fragment_options_str,
    ##                     index_options_str, int(index_options.symmetric))
        
    def add_dataset(self, title, creation_date, fragment_options_str,
                   index_options_str, is_symmetric):
        self.conn.execute(self.ADD_DATASET_SQL,
                          (title, creation_date, fragment_options_str,
                            index_options_str, is_symmetric))

    ## def _get_table_count(self, table):
    ##     assert table in schema.SCHEMA_TABLE_NAMES, table
    ##     sql = "SELECT COUNT(*) FROM %s" % (table,)
    ##     sql = self.prepare_single(sql)
    ##     self.conn.execute(sql)
    ##     return schema._get_one(self.conn)

    def flush(self):
        pass

    def end(self, reporter):
        self.flush()
        reporter.update("Building index ...")
        schema.create_index(self.conn)
        
        # Improve query planning
        reporter.update("Analyzing database ...")
        self.conn.execute("ANALYZE")
        
        reporter.update("Computing sizes ...")

        sql = self.prepare_single(UPDATE_DATASET_SQL)
        self.conn.execute(sql)
        
        sql = self.prepare_single(UPDATE_RULE_ENV_NUM_PAIRS)
        self.conn.execute(sql)
        
        reporter.update("")

    if 0:
        # Experimental method to auto-create a database
        @staticmethod
        def create_database(uri, connect_kwargs):
            pass
####

class SingleIndexWriterMixin(object):
    """Each add_*() call corresponds to a single database execute()"""
    def add_property_name(self, property_name_idx, property_name):
        self.conn.execute(self.ADD_PROPERTY_NAME_SQL,
                              (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self.conn.execute(self.ADD_RULE_SMILES_SQL,
                              (smiles_idx, smiles, get_num_heavies_from_smiles(smiles)))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self.conn.execute(self.ADD_RULE_SQL,
                              (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, smarts, pseudosmiles, parent_smarts):
        if parent_smarts is None:
            parent_smarts = ""
        self.conn.execute(self.ADD_ENVIRONMENT_FINGERPRINT_SQL,
                              (fp_idx, smarts, pseudosmiles, parent_smarts))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self.conn.execute(self.ADD_RULE_ENVIRONMENT_SQL,
                              (rule_env_idx, rule_idx, env_fp_idx, radius, 0))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self.conn.execute(self.ADD_COMPOUND_SQL,
                              (compound_idx, compound_id, input_smiles, normalized_smiles, num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self.conn.execute(self.ADD_CONSTANT_SMILES_SQL,
                              (smiles_idx, constant_smiles))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self.conn.execute(self.ADD_RULE_ENVIRONMENT_PAIR_SQL,
                              (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self.conn.execute(self.ADD_COMPOUND_PROPERTY_SQL,
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
        
        self.conn.execute(self.ADD_RULE_ENVIRONMENT_STATISTICS_SQL,
                              (rule_env_idx, property_name_idx, count, avg, std,
                                   kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value))


# Helper generator to yield a 1 once every n steps.
# Used to trigger a flush() in the BatchIndexWriterMixin.
def _trigger(n):
    assert n > 0
    x = range(n-1)
    while 1:
        for _ in x:
            yield 0
        yield 1

class BatchIndexWriterMixin(object):
    """The add_*() calls are batched up, and sent via execute_many()"""
    def __init__(self, uri, db, conn, title):
        super(BatchIndexWriterMixin, self).__init__(uri, db, conn, title)
        self._property_name_values = []
        self._rule_smiles_values = []
        self._rule_values = []
        self._environment_fingerprint_values = []
        self._rule_environment_values = []
        self._compound_values = []
        self._constant_smiles_values = []
        self._rule_environment_pair_values = []
        self._compound_property_values = []
        self._rule_environment_statistics_values = []
        self._check_flush = _trigger(100000)

    def execute_many(self, sql, values):
        self.conn.executemany(sql, values)
        
    def add_property_name(self, property_name_idx, property_name):
        self._property_name_values.append(
            (property_name_idx, property_name) )
        if next(self._check_flush):
            self.flush()
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self._rule_smiles_values.append(
            (smiles_idx, smiles, get_num_heavies_from_smiles(smiles)))
        if next(self._check_flush):
            self.flush()

                       
    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self._rule_values.append(
            (rule_idx, from_smiles_idx, to_smiles_idx))
        if next(self._check_flush):
            self.flush()

    def add_environment_fingerprint(self, fp_idx, smarts, pseudosmiles, parent_smarts):
        if parent_smarts is None:
            parent_smarts = ""
        self._environment_fingerprint_values.append(
                       (fp_idx, smarts, pseudosmiles, parent_smarts))
        if next(self._check_flush):
            self.flush()

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        # also initialize num_pairs to 0
        self._rule_environment_values.append(
            (rule_env_idx, rule_idx, env_fp_idx, radius, 0))
        if next(self._check_flush):
            self.flush()

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self._compound_values.append(
                       (compound_idx, compound_id, input_smiles,
                            normalized_smiles, num_normalized_heavies))
        if next(self._check_flush):
            self.flush()
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self._constant_smiles_values.append(
                       (smiles_idx, constant_smiles))
        if next(self._check_flush):
            self.flush()

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self._rule_environment_pair_values.append(
            (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))
        if next(self._check_flush):
            self.flush()
                 

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self._compound_property_values.append(
            (compound_idx, property_name_idx, value))
        if next(self._check_flush):
            self.flush()

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        self._rule_environment_statistics_values.append(
            (rule_env_idx, property_name_idx) + tuple(values))
        if next(self._check_flush):
            self.flush()

        
    def flush(self):
        for (values, sql) in (
                (self._property_name_values, self.ADD_PROPERTY_NAME_SQL),
                (self._rule_smiles_values, self.ADD_RULE_SMILES_SQL),
                (self._rule_values, self.ADD_RULE_SQL),
                (self._environment_fingerprint_values,
                     self.ADD_ENVIRONMENT_FINGERPRINT_SQL),
                (self._rule_environment_values, self.ADD_RULE_ENVIRONMENT_SQL),
                (self._compound_values, self.ADD_COMPOUND_SQL),
                (self._constant_smiles_values, self.ADD_CONSTANT_SMILES_SQL),
                (self._rule_environment_pair_values, self.ADD_RULE_ENVIRONMENT_PAIR_SQL),
                (self._compound_property_values,
                     self.ADD_COMPOUND_PROPERTY_SQL),
                (self._rule_environment_statistics_values,
                     self.ADD_RULE_ENVIRONMENT_STATISTICS_SQL)):
            if values:
                ## print("Execute many", repr(sql))
                ## for q in values:
                ##     print("GOT", repr(q))
                self.execute_many(sql, values)
                del values[:]
        

class NoTransactionMixin(object):
    """Don't wrap the indexing in BEGIN TRANSACTION / COMMIT"""
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
    

class TransactionMixin(object):
    """Wrap the indexing in BEGIN TRANSACTION / COMMIT"""
    def start(self, fragment_options, index_options):
        self.conn.execute("BEGIN TRANSACTION")
        super(TransactionMixin, self).start(fragment_options, index_options)
    
    def close(self):
        #self.db.execute("COMMIT")
        self.db.close()

    def commit(self):
        self.conn.execute("COMMIT")
        self.conn.close()
        self.db.close()

    def rollback(self):
        #self.conn.execute("ROLLBACK")
        self.conn.close()
        self.db.close()


class BaseSQLiteIndexWriter(SingleIndexWriterMixin, BaseRDBMSIndexWriter):
    """Base class for using Python's own 'sqlite3' module or the 'apsw' module"""
    def create_schema(self):
        schema.create_schema_for_sqlite(self.db)
        
class SQLiteIndexWriter(NoTransactionMixin, BaseSQLiteIndexWriter):
    pass

class APSWIndexWriter(TransactionMixin, BaseSQLiteIndexWriter):
    pass


## For Postgres support

def _get_main_postgres_database(uri):
    from playhouse import db_url

    try:
        from urlparse import urlparse, urlunparse
    except ImportError:
        from urllib.parse import urlparse, urlunparse

    parsed = urlparse(uri)
    new_parsed = parsed._replace(path="/postgres")
    postgres_uri = urlunparse(new_parsed)

    return db_url.connect(postgres_uri)

def _quote_postgres_database(conn, database_name):
    import psycopg2.sql
    return psycopg2.sql.Identifier(database_name).as_string(conn)


class PostgresIndexWriter(TransactionMixin, BatchIndexWriterMixin, BaseRDBMSIndexWriter):
    if 0:
        # Experimental method to auto-create a database
        @staticmethod
        def create_database(uri, connect_kwargs):
            # Create the database. Don't create the schema.
            import psycopg2
            import psycopg2.sql

            # If the database doesn't exist, then create it
            if "database" not in connect_kwargs:
                raise ValueError("Daabtase URI %r does not contain a database" % (uri,))
            database = connect_kwargs["database"]
            if database == "postgres":
                raise ValueError("Cannot use the Postgres database named 'postgres'")

            conn = _get_main_postgres_database(uri)
            conn.connect()
            c = conn.cursor()
            c.execute("SELECT datname FROM pg_database WHERE datname = %s and datistemplate = false;",
                          (database,))
            databases = [row[0] for row in c]
            if database not in databases:
                # Need to create the database. Peewee doesn't handle this so work
                # with the low-level psycopg2 interface.
                import psycopg2.sql

                raw_conn = conn.connection()

                # Need to change the isolation level in order to create the database
                level = raw_conn.isolation_level
                raw_conn.set_isolation_level(0)  # ISOLATION_LEVEL_AUTOCOMMIT
                try:

                    # Database names use a different format than row names,
                    # so can't use "%s" for that term.
                    sql = "CREATE DATABASE {} WITH ENCODING='UTF8' LC_COLLATE='C' LC_CTYPE='C'"
                    query = psycopg2.sql.SQL(sql.format(
                        _quote_postgres_database(raw_conn), connect_kwargs["database"]))
                    raw_conn.cursor().execute(query) # create the database
                finally:
                    raw_conn.set_isolation_level(level)

            c.close()
            conn.close()

    def create_schema(self):
        # Check that the database is empty
        c = self.conn
        c.execute("SELECT table_name FROM information_schema.tables")
        table_names = set(row[0] for row in c)
        
        existing_table_names = []
        for table_name in schema.SCHEMA_TABLE_NAMES:
            if table_name in table_names:
                existing_table_names.append(table_name)

        if existing_table_names:
            for table_name in existing_table_names:
                c.execute("DROP TABLE {} CASCADE".format(_quote_postgres_database(c, table_name)))
            else:
                from .index_algorithm import DatabaseAlreadyExists
                tables = ", ".join(existing_table_names)
                if len(existing_table_names) == 1:
                    msg = f"Required table already exists: {tables}"
                else:
                    msg = f"Requires tables already exist: {tables}"
                
                raise DatabaseAlreadyExists("Postgres database", self.uri, msg)
        
        schema.create_schema(self.db, schema.PostgresConfig)

    def prepare(self, sql):
        return sql.replace("?", "%s")

    # The SQL statements are specially constructed to make this post-processing simple.
    def prepare_many(self, sql):
        left, mid, right = sql.partition("VALUES")
        sql = left + mid + " %s"
        return sql
            
        
    # execute_values() is faster than executemany()
    def execute_many(self, sql, values):
        from psycopg2.extras import execute_values
        execute_values(self.conn, sql, values)


class FlatSQLFile(object):
    def __init__(self, outfile):
        self.outfile = outfile
        
    def cursor(self):
        return FlatSQLFileCursor(self.outfile)

    def close(self):
        self.outfile.close()
    
    def commit(self):
        self.close()

    def rollback(self):
        self.close()

    def execute(self, s):
        if s != "COMMIT":
            raise AssertionError("Unexpected execute: %r" % (s,))
        
        
_question_mark_pat = re.compile("[?]")
_end_with_semicolon = re.compile(r";\s*$")

_question_substitutions = dict(
    (i, ("(" + ", ".join("?"*i) + ")")) for i in range(30)
    )

class FlatSQLFileCursor(object):
    def __init__(self, outfile):
        self.outfile = outfile

    def execute(self, sql, args=None):
        if args:
            n = sql.count("?")
            assert n == len(args), (n, sql, args)
            substitution = _question_substitutions[n]
            new_args = []
            for arg in args:
                if arg is None:
                    new_args.append("NULL")
                elif isinstance(arg, (int, float)):
                    new_args.append(str(arg))
                else:
                    assert isinstance(arg, str), arg
                    new_args.append("'" + arg.replace("'", "''") + "'")
            assert n == len(new_args), (n, sql, args)
            if substitution not in sql:
                # Only happens with "INSERT INTO dataset ... VALUES (1, 2, ?, ?, ..."
                new_args.reverse()
                sql = _question_mark_pat.sub(lambda m: new_args.pop(), sql)
            else:
                sql = sql.replace(substitution, " ".join(new_args))
        if sql[-2:] == ";\n" or sql[-1:] == ";":
            self.outfile.write(sql)
        else:
            self.outfile.write(sql + ";\n")
                                    

    def executemany(self, sql, values):
        if not values:
            pass
        n = len(values[0])
        assert sql.count("?") == n, (n, sql)
        substitution = _question_substitutions[n]
        assert substitution in sql, (n, sql)
        
        new_values = []
        for args in values:
            assert len(args) == n, (n, sql, args)
            new_args = []
            for arg in args:
                if arg is None:
                    new_args.append("NULL")
                elif isinstance(arg, (int, float)):
                    new_args.append(str(arg))
                else:
                    assert isinstance(arg, str), arg
                    new_args.append("'" + arg.replace("'", "''") + "'")
            new_values.append("\n    (" + ", ".join(new_args) + ")")

        sql = sql.replace(substitution, ",".join(new_values))
        self.outfile.write(sql + ";\n")


    def close(self):
        pass

class BaseSQLWriter(TransactionMixin, BatchIndexWriterMixin, BaseRDBMSIndexWriter):
    pass

class SQLiteSQLWriter(BaseSQLWriter):
    def create_schema(self):
        schema.create_schema_for_sqlite(self.db)

class PostgresSQLWriter(BaseSQLWriter):
    def create_schema(self):
        schema.create_schema(self.db, schema.PostgresConfig)

    
def open_sql_writer(destination, title, variant, compression):
    outfile = _open_output(destination, compression)
    db = FlatSQLFile(outfile)
    if variant == "sqlite":
        writer_class = SQLiteSQLWriter
    elif variant == "postgres":
        writer_class = PostgresSQLWriter
    else:
        raise ValueError("variant must be one of 'sqlite' or 'postgres'")
        
    conn = db.cursor()
    writer = writer_class(destination, db, conn, title)
    writer.create_schema()
    return writer

####

def open_rdbms_index_writer(filename, title, is_sqlite=False):
    # See if this is a database URI
    from playhouse import db_url
    import peewee

    if is_sqlite:
        database_class = None
    else:
        parsed = db_url.urlparse(filename)
        if parsed.scheme == "oracle":
                # peewee doesn't even support Oracle
                raise NotImplementedError("Oracle not supported")

        database_class = db_url.schemes.get(parsed.scheme, None)
    
    if database_class is not None:
        connect_kwargs = db_url.parseresult_to_dict(parsed)
        
        if issubclass(database_class, peewee.SqliteDatabase):
            writer_class = SQLiteIndexWriter
        elif issubclass(database_class, peewee.PostgresqlDatabase):
            writer_class = PostgresIndexWriter
        elif issubclass(database_class, peewee.MySQLDatabase):
            # I haven't developed this
            raise NotImplementedError("MySQL not supported")
        else:
            raise NotImplementedError("%r (%r) not supported" % (parsed.scheme, database_class))

        ## Experimental: Create database if it doesn't exist.
        #writer_class.create_database(filename, connect_kwargs)
        db = database_class(**connect_kwargs)
        db.connect()

    else:
        # Filename. Use the SQLite interface.

        if filename != ":memory:":
            if os.path.exists(filename):
                os.unlink(filename)
        if apsw is None:
            db = sqlite3.connect(filename)
            writer_class = SQLiteIndexWriter
        else:
            db = apsw.Connection(filename)
            writer_class = APSWIndexWriter
    
    conn = db.cursor()
    
    writer = writer_class(filename, db, conn, title)
    writer.create_schema()
    return writer

def update_counts(cursor):
    cursor.execute(UPDATE_DATASET_SQL)

def update_env_fp_num_pairs(cursor):
    cursor.execute(UPDATE_ENV_FP_NUM_PAIRS)
    
### format-specific writers
def _open_output(destination, compression):
    if compression:
        format_hint = ".gz"
    else:
        format_hint = ""
    return fileio.open_output(destination, format_hint)
    

# Can be a helpful summary
def _open_csv(destination, compression,
                title, fragment_options, fragment_index, index_options, properties,
                environment_cache):
    outfile = _open_output(destination, compression)
    return CSVPairWriter(outfile, fragment_options, fragment_index,
                         index_options, properties)

# Text-based version somewhat useful for debugging
def _open_mmpa(destination, compression, title):
    outfile = _open_output(destination, compression)
    return open_table_index_writer(outfile)

# SQLite database
def _open_mmpdb(destination, compression, title):
    if compression:
        raise ValueError("mmpdb output does not support compression")
    return open_rdbms_index_writer(destination, title)

# SQL for SQLite 
def _open_sqlite_sql(destination, compression, title):
    return open_sql_writer(destination, title, "sqlite", compression)

# SQL for Postgres
def _open_postgres_sql(destination, compression, title):
    return open_sql_writer(destination, title, "postgres", compression)

def NO_TABS(s):
    return s.replace("\t", " ")

def NULLABLE(s):
    if s is None:
        return "\\N" # Postgres syntax
    return str(s)

def write_readme(dirname):
    with open(os.path.join(dirname, "README"), "w") as f:
        f.write("""\
CSV files for an mmpdb matched molecular pair database, along with
drivers to import the data into SQLite and Postgres.

To import the dataset into the sqlite database 'DB.mmpdb' use the
following:

   sqlite3 DB.mmpdb < load_mmpdb.sqlite

The sqlite loader 'load_mmpdb.sqlite':
  1) creates the schema using 'load_mmpdb_schema.sqlite'
  2) imports the data using 'load_mmpdb_data.sqlite'
  3) create indices using 'create_index.sql'
  4) analyzes the database

To import the dataset into the postgres database 'DB' use the
following:

   pgsql DB < load_mmpdb.psql

The sqlite loader 'load_mmpdb.sqlite':
  1) creates the schema using 'load_mmpdb_schema.psql'
  2) imports the data using 'load_mmpdb_data.psql'
  3) create indices using 'create_index.sql'
  4) analyzes the database
""")

def write_index_files(dirname):
    create_index_sql = schema.get_create_index_sql()
    filename = os.path.join(dirname, "create_index.sql")
    with open(filename, "w") as f:
        f.write(create_index_sql)


LOAD_SQLITE_SCRIPT = """\
.mode csv
.separator "\t"
.import compound.csv compound
.import compound_property.csv compound_property
.import constant_smiles.csv constant_smiles
.import dataset.csv dataset
.import environment_fingerprint.csv environment_fingerprint
.import property_name.csv property_name
.import rule.csv rule
.import rule_environment.csv rule_environment
.import pair.csv pair
.import rule_environment_statistics.csv rule_environment_statistics
.import rule_smiles.csv rule_smiles

"""

SQLITE_FIX_NULLS = """\
UPDATE rule_environment_statistics SET avg = NULL WHERE avg = '\\N';
UPDATE rule_environment_statistics SET std = NULL WHERE std = '\\N';
UPDATE rule_environment_statistics SET kurtosis = NULL WHERE kurtosis = '\\N';
UPDATE rule_environment_statistics SET skewness = NULL WHERE skewness = '\\N';
UPDATE rule_environment_statistics SET paired_t = NULL WHERE paired_t = '\\N';
UPDATE rule_environment_statistics SET p_value = NULL WHERE p_value = '\\N';
"""

def write_sqlite_files(dirname):
    with open(os.path.join(dirname, "load_mmpdb.sqlite"), "w") as f:
        f.write(".read load_mmpdb_schema.sqlite\n")
        f.write(".read load_mmpdb_data.sqlite\n")
        f.write(".read create_index.sql\n")
        f.write("ANALYZE\n")

    filename = os.path.join(dirname, "load_mmpdb_schema.sqlite")
    open_sql_writer(filename, None, "sqlite", False).close()
        
    with open(os.path.join(dirname, "load_mmpdb_data.sqlite"), "w") as f:
        f.write(LOAD_SQLITE_SCRIPT)
        f.write(UPDATE_DATASET_SQL + ";\n")
        f.write(UPDATE_RULE_ENV_NUM_PAIRS + ";\n")
        f.write(SQLITE_FIX_NULLS)

LOAD_PSQL_SCRIPT = r"""
\echo Loading 'dataset' table...
-- The 'sed' commands are used to escape backslashes appropriate for Postgres.
-- They are only needed for tables which may contain strings with a backslash.
\copy dataset FROM PROGRAM 'sed ''s/\\[\N]/\\\\/g'' < dataset.csv' WITH (FORMAT text);

\echo Loading 'property_name' table...
\copy property_name FROM PROGRAM 'sed ''s/\\/\\\\/g'' < property_name.csv' WITH (FORMAT text);

\echo Loading 'rule_smiles' table...
\copy rule_smiles FROM PROGRAM 'sed ''s/\\/\\\\/g'' < rule_smiles.csv' WITH (FORMAT text);

\echo Loading 'rule' table...
\copy rule FROM 'rule.csv' WITH (FORMAT text);

\echo Loading 'environment_fingerprint' table...
\copy environment_fingerprint FROM 'environment_fingerprint.csv' WITH (FORMAT text);

\echo Loading 'rule_environment' table...
\copy rule_environment FROM 'rule_environment.csv' WITH (FORMAT text);

\echo Loading 'compound' table...
\copy compound FROM PROGRAM 'sed ''s/\\/\\\\/g'' < compound.csv' WITH (FORMAT 'text');

\echo Loading 'constant_smiles' table...
\copy constant_smiles FROM PROGRAM 'sed ''s/\\/\\\\/g'' < constant_smiles.csv' WITH (FORMAT 'text');

\echo Loading 'pair' table...
\copy pair FROM 'pair.csv' WITH (FORMAT text);

\echo Loading 'compound_property' table...
\copy compound_property FROM 'compound_property.csv' WITH (FORMAT 'text');

\echo Loading 'rule_environment_statistics' table...
-- Note: some of the fields may be \N which must be interpreted as a NULL.
-- Do not escape it with the simple sed command. (No need, as there are no strings.)
\copy rule_environment_statistics FROM 'rule_environment_statistics.csv' WITH (FORMAT text);

\echo Updating counts ...
"""
        
def write_postgres_files(dirname):
    with open(os.path.join(dirname, "load_mmpdb.psql"), "w") as f:
        f.write("\\echo Creating mmpdb schema.\n")
        f.write("\\i load_mmpdb_schema.psql\n")
        f.write("\\echo Copying mmpdb dataset to postgres server.\n")
        f.write("\\i load_mmpdb_data.psql\n")
        f.write("\\echo Creating index.\n")
        f.write("\\i create_index.sql\n")
        f.write("\\echo Analyzing.\n")
        f.write("ANALYZE\n")
        f.write("\\echo Done.\n")

    filename = os.path.join(dirname, "load_mmpdb_schema.psql")
    open_sql_writer(filename, None, "postgres", False).close()
        
    with open(os.path.join(dirname, "load_mmpdb_data.psql"), "w") as f:
        f.write(LOAD_PSQL_SCRIPT)
        f.write(UPDATE_DATASET_SQL + ";\n")
        f.write(UPDATE_RULE_ENV_NUM_PAIRS + ";\n")
        f.write("\\echo Finished copying mmpdb dataset.\n")


class CSVDWriter(StartMixin):
    def __init__(self,
                     destination,
                     title,
                     dataset_file,
                     compound_file,
                     property_name_file,
                     compound_property_file,
                     rule_smiles_file,
                     rule_file,
                     environment_fingerprint_file,
                     rule_environment_file,
                     pair_file,
                     constant_smiles_file,
                     rule_environment_statistics_file):
        self.destination = destination
        self.title = title
        
        self.dataset_file = dataset_file
        self.compound_file = compound_file
        self.property_name_file = property_name_file
        self.compound_property_file = compound_property_file
        self.rule_smiles_file = rule_smiles_file
        self.rule_file = rule_file
        self.environment_fingerprint_file = environment_fingerprint_file
        self.rule_environment_file = rule_environment_file
        self.pair_file = pair_file
        self.constant_smiles_file = constant_smiles_file
        self.rule_environment_statistics_file = rule_environment_statistics_file

        self._rule_environment_statistics_id_gen = itertools.count(1)
        self._compound_property_id_gen = itertools.count(1)

    def start(self, *args, **kwargs):
        super(CSVDWriter, self).start(*args, **kwargs)
        write_readme(self.destination)
        write_index_files(self.destination)
        write_sqlite_files(self.destination)
        write_postgres_files(self.destination)
        

    def rollback(self):
        # XXX should I delete the files?
        self.close()

    def end(self, reporter):
        pass

    def commit(self):
        self.close()
    
    def close(self):
        self.dataset_file.close()
        self.compound_file.close()
        self.property_name_file.close()
        self.compound_property_file.close()
        self.rule_smiles_file.close()
        self.rule_file.close()
        self.environment_fingerprint_file.close()
        self.rule_environment_file.close()
        self.pair_file.close()
        self.constant_smiles_file.close()
        self.rule_environment_statistics_file.close()
    
    def add_dataset(self, title, creation_date, fragment_options_str,
                   index_options_str, is_symmetric):
        # Fill in num_compounds, num_rules, ... num_rule_environment_stats with -1
        self.dataset_file.write("1\t4\t%s\t%s\t%s\t%s\t%d\t-1\t-1\t-1\t-1\t-1\n" % (
            (title, NO_TABS(creation_date), NO_TABS(fragment_options_str),
                 NO_TABS(index_options_str), is_symmetric)))
        
    def add_property_name(self, property_name_idx, property_name):
        self.property_name_file.write("%d\t%s\n" % (property_name_idx, NO_TABS(property_name)))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self.rule_smiles_file.write("%d\t%s\t%d\n" % (smiles_idx, NO_TABS(smiles),
                                                  get_num_heavies_from_smiles(smiles)))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self.rule_file.write("%d\t%d\t%d\n" % (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, smarts, pseudosmiles, parent_smarts):
        if parent_smarts is None:
            parent_smarts = ""
        self.environment_fingerprint_file.write("%d\t%s\t%s\t%s\n" % (
            fp_idx, smarts, pseudosmiles, parent_smarts))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self.rule_environment_file.write("%d\t%d\t%d\t%d\n" % (rule_env_idx, rule_idx, env_fp_idx, radius))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self.compound_file.write("%d\t%s\t%s\t%s\t%d\n" % 
                 (compound_idx, NO_TABS(compound_id), NO_TABS(input_smiles),
                      NO_TABS(normalized_smiles), num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self.constant_smiles_file.write("%d\t%s\n" % (smiles_idx, NO_TABS(constant_smiles)))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self.pair_file.write("%d\t%d\t%d\t%d\t%d\n" %
                              (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self.compound_property_file.write("%d\t%d\t%d\t%s\n" %
                              (next(self._compound_property_id_gen),
                                   compound_idx, property_name_idx, value))

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = values
        self.rule_environment_statistics_file.write(
            "%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
            (next(self._rule_environment_statistics_id_gen),
             rule_env_idx, property_name_idx, count, avg, NULLABLE(std),
             NULLABLE(kurtosis), NULLABLE(skewness), min, q1, median, q3, max,
             NULLABLE(paired_t), NULLABLE(p_value)))


def _open_csvd(destination, compression, title):
    os.makedirs(destination, exist_ok=True)
    writer = CSVDWriter(
        destination = destination,
        title = title,
        dataset_file = open(os.path.join(destination, "dataset.csv"), "w"),
        compound_file = open(os.path.join(destination, "compound.csv"), "w"),
        property_name_file = open(os.path.join(destination, "property_name.csv"), "w"),
        compound_property_file = open(os.path.join(destination, "compound_property.csv"), "w"),
        rule_smiles_file = open(os.path.join(destination, "rule_smiles.csv"), "w"),
        rule_file = open(os.path.join(destination, "rule.csv"), "w"),
        environment_fingerprint_file = open(os.path.join(destination, "environment_fingerprint.csv"), "w"),
        rule_environment_file = open(os.path.join(destination, "rule_environment.csv"), "w"),
        pair_file = open(os.path.join(destination, "pair.csv"), "w"),
        constant_smiles_file = open(os.path.join(destination, "constant_smiles.csv"), "w"),
        rule_environment_statistics_file = open(os.path.join(destination, "rule_environment_statistics.csv"), "w"),
        )
    return writer
    

class BaseFormatInfo(object):
    def __init__(self, format_name, supports_compression, opener, help, extensions=None):
        self.format_name = format_name
        if extensions is None:
            extensions = [format_name]
            if supports_compression:
                extensions.append(format_name + ".gz")
        self.supports_compression = supports_compression
        self.opener = opener
        self.help = help
        self.extensions = extensions

    def open(self, destination, compression,
                 title, fragment_options, fragment_index, index_options, properties,
                 environment_cache):
        return self.opener(destination, compression,
                               title, fragment_options, fragment_index, index_options, properties,
                               environment_cache)
            
class FormatInfo(BaseFormatInfo):
    def open(self, destination, compression,
                 title, fragment_options, fragment_index, index_options, properties,
                 environment_cache):
        index_writer = self.opener(destination, compression, title)
        return index_algorithm.MMPWriter(index_writer, fragment_options, fragment_index,
                                             index_options, properties)
        
# format name, [list, of, extensions]
_format_info_list = [
    FormatInfo("mmpdb", False, _open_mmpdb, "SQLite database"),
    BaseFormatInfo("csv", True, _open_csv, "simple CSV output"),
    FormatInfo("mmpa", True, _open_mmpa, "simple text format"),
    FormatInfo("sqlite", True, _open_sqlite_sql, "SQL schema and data dump, for SQLite"),
    FormatInfo("sql", True, _open_sqlite_sql, "alias for 'sqlite'"),
    FormatInfo("postgres", True, _open_postgres_sql, "SQL schema and data dump, for Postgres"),
    FormatInfo("csvd", False, _open_csvd, "mmpdb tables as CSV files and SQL scripts"),
]

_format_info_table = dict()
_extension_table = dict()

def _init_format_info():
    for format_info in _format_info_list:
        _format_info_table[format_info.format_name] = format_info
        
        for ext in format_info.extensions:
            _extension_table["." + ext] = format_info

_init_format_info()

def open_mmpa_writer(destination, format, title, fragment_options,
                     fragment_index, index_options, properties,
                     environment_cache):
    compression = False
    
    if format is None:
        if destination is None:
            format_info = _format_info_table["mmpa"]
        else:
            filename = destination.lower()
            if filename.endswith(".gz"):
                compression = True
                filename = filename[:-3]
            ext = os.path.splitext(filename)[1]
            format_info = _extension_table.get(ext, None)

            if format_info is None:
                format_info = _format_info_table["mmpdb"]
                
            if compression and not format_info.supports_compression:
                raise ValueError("Compression not available for %r" % (destination,))

    else:
        # Format specified. See if it's supported.
        input_format = format
        if format.endswith(".gz"):
            compression = True
            format = format[:-3]
        format_info = _format_info_table.get(format, None)
        if format_info is None:
            raise ValueError("Unsupported format: %r" % (input_format,))
        
        if compression and not format_info.supports_compression:
                raise ValueError("Compression not available for format %r" % (input_format,))

    return format_info.open(destination, compression, title, fragment_options,
                                fragment_index, index_options, properties,
                                environment_cache)
