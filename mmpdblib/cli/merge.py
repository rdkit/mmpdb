import os
import time
import sqlite3

import click

from .click_utils import (
    command,
    die,
    add_multiple_databases_parameters,
    )

from .. import dbutils
from .. import index_writers
from .. import schema


LUT_CREATE_TABLE = """

CREATE TABLE lut.rule_map (
  old_id INTEGER PRIMARY KEY,
  new_id INTEGER
);

CREATE TABLE lut.env_fp_map (
  old_id INTEGER PRIMARY KEY,
  new_id INTEGER
);

CREATE TABLE lut.rule_env_map (
  old_id INTEGER PRIMARY KEY,
  new_id INTEGER
);

"""

MERGE_CREATE_INDEX_SQL = """
CREATE UNIQUE INDEX idx_compound ON compound(public_id);
CREATE UNIQUE INDEX idx_rule_smiles ON rule_smiles(smiles);
CREATE UNIQUE INDEX idx_rule ON rule (from_smiles_id, to_smiles_id);
CREATE UNIQUE INDEX idx_env_fp ON environment_fingerprint (smarts, parent_smarts);
CREATE UNIQUE INDEX idx_rule_env ON rule_environment (rule_id, environment_fingerprint_id, radius);
CREATE UNIQUE INDEX idx_constant_smiles ON constant_smiles(smiles);
"""

MERGE_DROP_INDEX_SQL = """
DROP INDEX idx_compound;
DROP INDEX idx_rule_smiles;
DROP INDEX idx_rule;
DROP INDEX idx_env_fp;
DROP INDEX idx_rule_env;
DROP INDEX idx_constant_smiles;
"""

# Here's the schema to create the --profile database
PROFILE_SCHEMA = """
CREATE TABLE IF NOT EXISTS profiler.mmpdb_times (
  id integer primary key,
  label string,
  t timestamp
);
"""

MERGE_DATABASE_SQL = """

-- This expects the database to import to be attached as 'old'
-- using something like:
--   ATTACH DATABASE "subset.000.mmpdb" AS old
-- 
-- and a timer database (which may be ":memory:") as 'profiler', like:
--   ATTACH DATABASE "times.db" as profiler;

-- Step 1: Copy over the compound table

-- Must have UNIQUE INDEX ON compound (public_id)
-- The merge code ensured duplicate entries have the same clean_smiles.


INSERT INTO profiler.mmpdb_times (label, t) VALUES ("start merge", julianday("now"));

INSERT OR IGNORE INTO compound (public_id, input_smiles, clean_smiles, clean_num_heavies)
 SELECT public_id, input_smiles, clean_smiles, clean_num_heavies
   FROM old.compound
        ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert compound", julianday("now"));

-- Step 2: Copy over the rule_smiles

-- Must have UNIQUE INDEX ON rule_smiles (smiles)

INSERT OR IGNORE INTO rule_smiles (smiles, num_heavies)
 SELECT smiles, num_heavies
   FROM old.rule_smiles
        ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert rule_smiles", julianday("now"));

-- Step 3: Merge any rules

-- Must have UNIQUE INDEX ON rule (from_smiles_id, to_smiles_id)

INSERT OR IGNORE INTO rule (from_smiles_id, to_smiles_id)
 SELECT new_from_smiles.id, new_to_smiles.id
   FROM old.rule AS old_rule,
        old.rule_smiles AS old_from_smiles,
        old.rule_smiles AS old_to_smiles,
        rule_smiles AS new_from_smiles,
        rule_smiles AS new_to_smiles
  WHERE (    old_rule.from_smiles_id = old_from_smiles.id
         AND old_rule.to_smiles_id = old_to_smiles.id)
    AND (old_from_smiles.smiles = new_from_smiles.smiles
         AND old_to_smiles.smiles = new_to_smiles.smiles)
        ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert rule", julianday("now"));

-- Step 4. Table mapping old rule id to new rule id

-- This simplifies steps 6 and 7

DELETE FROM lut.rule_map;

INSERT INTO lut.rule_map (old_id, new_id) 
 SELECT old_rule.id, new_rule.id
   FROM old.rule AS old_rule,
        old.rule_smiles AS old_from_smiles,
        old.rule_smiles AS old_to_smiles,
        rule AS new_rule,
        rule_smiles AS new_from_smiles,
        rule_smiles AS new_to_smiles
  WHERE (old_rule.from_smiles_id = old_from_smiles.id
         AND old_rule.to_smiles_id = old_to_smiles.id)
    AND (old_from_smiles.smiles = new_from_smiles.smiles
         AND old_to_smiles.smiles = new_to_smiles.smiles)
    AND (new_rule.from_smiles_id = new_from_smiles.id
         AND new_rule.to_smiles_id = new_to_smiles.id)
        ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert rule_map", julianday("now"));

-- Step 5: Merge environment_fingerprint

-- Must have UNIQUE INDEX ON environment_fingerprint (smarts, parent_smarts)
-- (I haven't checked if I only need to use the smarts.)

INSERT OR IGNORE INTO environment_fingerprint (smarts, pseudosmiles, parent_smarts)
  SELECT smarts, pseudosmiles, parent_smarts
    FROM old.environment_fingerprint
         ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert env_fp", julianday("now"));

DELETE FROM lut.env_fp_map;

INSERT INTO lut.env_fp_map (old_id, new_id)
  SELECT old_env_fp.id, new_env_fp.id
    FROM old.environment_fingerprint as old_env_fp,
         environment_fingerprint as new_env_fp
   WHERE old_env_fp.smarts = new_env_fp.smarts
     AND old_env_fp.parent_smarts = new_env_fp.parent_smarts
         ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert env_fp_map", julianday("now"));

-- Step 6: Merge rule_environment

-- Must have UNIQUE INDEX ON rule_environment (rule_id, environment_fingerprint_id, radius)


INSERT OR IGNORE INTO rule_environment (rule_id, environment_fingerprint_id, radius)
  SELECT rule_map.new_id, env_fp_map.new_id, old_rule_env.radius
    FROM old.rule_environment AS old_rule_env,
         lut.rule_map AS rule_map,
         lut.env_fp_map AS env_fp_map
   WHERE old_rule_env.rule_id = rule_map.old_id
     AND old_rule_env.environment_fingerprint_id = env_fp_map.old_id
         ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert rule_environment", julianday("now"));

-- Step 7: Table mapping old rule environment to new rule environment


DELETE FROM lut.rule_env_map;

INSERT INTO lut.rule_env_map (old_id, new_id)
  SELECT old_rule_env.id, new_rule_env.id
    FROM old.rule_environment AS old_rule_env,
         rule_environment AS new_rule_env,
         lut.rule_map AS rule_map,
         lut.env_fp_map AS env_fp_map
   WHERE old_rule_env.rule_id = rule_map.old_id
     AND new_rule_env.rule_id = rule_map.new_id
     AND old_rule_env.environment_fingerprint_id = env_fp_map.old_id
     AND new_rule_env.environment_fingerprint_id = env_fp_map.new_id
     AND old_rule_env.radius = new_rule_env.radius
         ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert rule_env_map", julianday("now"));

-- Step 8: Merge constant_smiles

-- requires a unique index

INSERT OR IGNORE INTO constant_smiles (smiles)
   SELECT smiles
     FROM old.constant_smiles
          ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert constant_smiles", julianday("now"));

-- Step 9: Merge pairs

-- These are all unique because the constant_smiles wasn't seen earlier

INSERT INTO pair (rule_environment_id, compound1_id, compound2_id, constant_id) 
 SELECT rule_env_map.new_id, new_compound1.id, new_compound2.id, new_constant.id
   FROM old.pair AS old_pair,
        lut.env_fp_map AS rule_env_map,
        old.compound as old_compound1,
        old.compound as old_compound2,
        compound as new_compound1,
        compound as new_compound2,
        old.constant_smiles as old_constant,
        constant_smiles as new_constant
  WHERE (    old_pair.rule_environment_id = rule_env_map.old_id
         AND old_pair.compound1_id = old_compound1.id
         AND old_pair.compound2_id = old_compound2.id
         AND old_pair.constant_id = old_constant.id)
    AND old_compound1.public_id = new_compound1.public_id
    AND old_compound2.public_id = new_compound2.public_id
    AND old_constant.smiles = new_constant.smiles
        ;

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("insert pairs", julianday("now"));

INSERT INTO profiler.mmpdb_times (label, t) VALUES ("end merge", julianday("now"));

"""

def _check_options_mismatch(category, dbinfo, options, first_dbinfo, first_options):
    d = options.to_dict()
    first_d = first_options.to_dict()
    if d == first_d:
        return

    # Figure out which values are different
    lines = [
        f"Cannot merge. The fragment options in {dbinfo.get_human_name()} "
        f"differ from {first_dbinfo.get_human_name()!r}."
        ]
    for k in d:
        if d[k] != first_d[k]:
            lines.append(f"  {k}: {d[k]!r} != {first_d[k]!r}")
    die(*lines)


def check_options_mismatch(dbinfo, fragment_options, index_options,
                           first_dbinfo, first_fragment_options, first_index_options):
    _check_options_mismatch("fragment", dbinfo, fragment_options, first_dbinfo, first_fragment_options)
    _check_options_mismatch("index", dbinfo, fragment_options, first_dbinfo, first_fragment_options)

def open_output_database(
        output_filename,
        title,
        fragment_options,
        index_options,
        replace,
        ):
    import sqlite3
    from .. import index_writers, reporters
    try:
        writer = index_writers.open_rdbms_index_writer(
            output_filename, title, replace=True, is_sqlite=True)
    except IOError as err:
        die(f"Cannot open SQLite file {output_filename!r}: {err}")
    except sqlite3.OperationalError as err:
        die(f"Unexpected problem initializing mmpdb file {output_filename!r}: {err}")

    writer.start(fragment_options, index_options)

    # Create the index
    writer.end(reporters.Quiet())
    writer.close()
    try:
        db = sqlite3.connect(output_filename)
    except sqlite3.OperationalError as err:
        die(f"Unexpected problem re-opening mmpdb file {output_filename!r}: {err}")

    #db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA synchronous=OFF")
    return db

def add_time(output_c, label):
    output_c.execute(
        'INSERT INTO profiler.mmpdb_times (label, t) VALUES (?, julianday("now"))',
        (label,),
        )

def attach_lut(output_c):
    output_c.execute('ATTACH DATABASE ":memory:" AS lut')
    schema._execute_sql(output_c, LUT_CREATE_TABLE)
    
def attach_profiler(output_c, profile_path):
    if profile_path is None:
        output_c.execute('ATTACH DATABASE ":memory:" AS profiler')
    else:
        try:
            output_c.execute('ATTACH DATABASE ? AS profiler', (profile_path,))
        except sqlite3.DatabaseError as err:
            die(f"Cannot attach profiler database {profile_path!r}: {err}")

    try:
        output_c.execute(PROFILE_SCHEMA)
    except sqlite3.OperationalError as err:
        die(f"Cannot create profiler schema in {profile_path!r}: {err}")
    
    add_time(output_c, "start")
    

def check_for_duplicate_constants(output_c, database):
    for smiles, in output_c.execute("""
SELECT old_constant_smiles.smiles
  FROM old.constant_smiles as old_constant_smiles,
       constant_smiles as new_constant_smiles
 WHERE old_constant_smiles.smiles = new_constant_smiles.smiles
"""):
        die(
            "Can only merge mmpdb databases with no shared constants.",
            f"Duplicate constant {smiles!r} found in {database!r} ."
            )


def check_for_different_clean_smiles(output_c, database):
    # Check any records with a different clean_smiles
    for public_id, old_clean_smiles, new_clean_smiles in output_c.execute("""
SELECT old_compound.public_id, old_compound.clean_smiles, new_compound.clean_smiles
  FROM old.compound as old_compound,
       compound as new_compound
 WHERE old_compound.public_id = new_compound.public_id
   AND old_compound.clean_smiles != new_compound.clean_smiles
"""):
        die(
            "Cannot merge. Record {public_id!r} in {database!r} has a SMILES {old_clean_smiles!r}",
            "but a previously merged file has the SMILES {new_clean_smiles!r}.",
            )

        
@command()
@click.option(
    "--title",
    help = "title to use for the output file",
    )
@click.option(
    "--output",
    "-o",
    "output_filename",
    default = None,
    type = click.Path(),
    )
@click.option(
    "--replace / --no-replace",
    default = False,
    help = "With --replace, replace any existing database. Default is --no-replace.",
    )
@click.option(
    "--profile",
    "profile_path",
    type = click.Path(),
    default = None,
    help = "store SQL profile information to this SQLite database",
    )
@add_multiple_databases_parameters()
@click.pass_obj
def merge(
        reporter,
        databases_options,
        title,
        output_filename,
        replace,
        profile_path,
        ):
    """merge multiple mmpdb databases

    Each DATABASE must be an mmpdb file.

    Properties are ignored because merging them is not meaningful.
    """
    
    databases = databases_options.databases
    if not databases:
        die("Must specify at least one mmpdb database.")

    if title is None:
        title = f"Merged MMPs from {len(databases)} files"

    if output_filename is None:
        output_filename = "merged.mmpdb"
        reporter.warning(f"No --output file specified, using {output_filename!r}")
    
    first_dbinfo = None
    first_fragment_options = None
    first_index_options = None

    output_db = None
    output_c = None

    db = None
    c = None
    try:
        for database_i, database in enumerate(databases, 1):
            # Require it to be a SQLite database
            if not os.path.exists(database):
                die(f"mmpdb file {database!r} does not exist")
            dbinfo = dbutils.DBFile(database)
            
            try:
                dataset = dbinfo.open_database(
                    quiet=reporter.quiet,
                    apsw_warning=False,
                    ).get_dataset()
            except dbutils.DBError as err:
                die(f"Cannot connect to {dbinfo.get_human_name()}: {err}\n")

            db = dataset.mmpa_db
            with db:
                # Get options
                fragment_options = dataset.get_fragment_options()
                index_options = dataset.get_index_options()

                # Warn that properties will not be merged
                property_names = dataset.get_property_names()
                if property_names:
                    msg = ", ".join(repr(name) for name in property_names)
                    reporter.warning(
                        f"Unable to merge properties from {dbinfo.get_human_name()}: {msg}"
                        )

            
                # Ensure the options are compatible
                if first_dbinfo is None:
                    # First time through
                    first_dbinfo = dbinfo
                    first_fragment_options = fragment_options
                    first_index_options = index_options

                    # Create the output database
                    output_db = open_output_database(
                        output_filename,
                        title,
                        fragment_options,
                        index_options,
                        replace,
                        )
                    output_c = output_db.cursor()
                    try:
                        schema._execute_sql(output_c, MERGE_CREATE_INDEX_SQL)
                    except sqlite3.DatabaseError as err:
                        output_c.close()
                        output_db.close()
                        output_c = output_db = None
                        die(f"Cannot use {output_filename!r}: {err}")
                        
                    except Exception:
                        output_c.close()
                        output_db.close()
                        output_c = output_db = None
                        raise

                    attach_lut(output_c)
                    attach_profiler(output_c, profile_path)

                else:
                    check_options_mismatch(
                            dbinfo,
                            fragment_options,
                            index_options,
                            first_dbinfo,
                            first_fragment_options,
                            first_index_options,
                            )

            # Attach and merge
            add_time(output_c, f"merging {database!r}")
            start_time = time.time()
            reporter.update(f"Merging {database!r} ({database_i}/{len(databases)}) ...")
            
            try:
                output_c.execute("ATTACH DATABASE ? AS old", (database,))
            except sqlite3.OperationalError as err:
                reporter.report(f"Merging {database!r}. FAILED!")
                die(f"Cannot attach {database!r} to {output_file!r}: {err}")

            try:
                try:
                    check_for_duplicate_constants(output_c, database)
                    check_for_different_clean_smiles(output_c, database)

                    # Merge
                    try:
                        schema._execute_sql(output_c, MERGE_DATABASE_SQL)
                    except Exception:
                        ## breakpoint()
                        raise
                    ## Debug code to see the query plan
##                     output_c.execute("""
## EXPLAIN
## INSERT OR IGNORE INTO rule_environment (rule_id, environment_fingerprint_id, radius)
##   SELECT rule_map.new_id, env_fp_map.new_id, old_rule_env.radius
##     FROM old.rule_environment AS old_rule_env,
##          lut.rule_map AS rule_map,
##          lut.env_fp_map AS env_fp_map
##    WHERE old_rule_env.rule_id = rule_map.old_id
##      AND old_rule_env.environment_fingerprint_id = env_fp_map.old_id
##          ;
## """)
##                     def STR(s, n):
##                         if s is None:
##                             s = ""
##                         return str(s).ljust(n)
##                     print()
##                     print("addr  opcode         p1    p2    p3    p4             p5  comment")
##                     print("----  -------------  ----  ----  ----  -------------  --  -------------")
##                     for addr, opcode, p1, p2, p3, p4, p5, comment in output_c:
##                         print(
##                             f"{STR(addr, 4)}  {STR(opcode, 13)}  {STR(p1, 4)}  {STR(p2, 4)}  "
##                             f"{STR(p3, 4)}  {STR(p4, 13)}  {STR(p5, 2)}  {STR(comment, 15)}"
##                             )
##                     print()
                            
                finally:
                    try:
                        output_c.execute("COMMIT")
                    except sqlite3.OperationalError:
                        pass
                    output_c.execute("DETACH DATABASE old")
                    # Following the recommendation at
                    #   https://sqlite.org/lang_analyze.html
                    output_c.execute("PRAGMA analysis_limit=400")
                    output_c.execute("PRAGMA optimize")
            except:
                reporter.report(f"Merging {database!r}. FAILED!")
                output_c.close()
                output_db.close()
                output_c = output_db = None
                raise
            reporter.report(
                f"Merged {database!r} ({database_i}/{len(databases)}). "
                f"Time: {time.time()-start_time:.2f}"
                )
                

    finally:
        if output_c is not None:
            try:
                output_c.execute("COMMIT")
            except sqlite3.OperationalError:
                pass
            output_c.execute("DETACH DATABASE lut")

            add_time(output_c, "start finalization")
            schema._execute_sql(output_c, MERGE_DROP_INDEX_SQL)
            add_time(output_c, "drop indices")

            index_writers.update_counts(output_c)
            add_time(output_c, "update counts")
            output_c.execute("ANALYZE")
            add_time(output_c, "analyze")
            add_time(output_c, "end finalization")
            add_time(output_c, "end")
            output_c.execute("COMMIT")
            output_c.close()
            output_db.close()

@command(
    name="show_merge_profile",
    )
@click.argument(
    "profile_path",
    type = click.Path(),
    metavar = "PROF.db",
    )
def show_merge_profile(profile_path):
    "Show information from a merge --profile database"
    if not os.path.exists(profile_path):
        die(f"Profile database does not exist: {profile_path!r}")

    try:
        db = sqlite3.connect(profile_path)
    except sqlite3.DatabaseError as err:
        die(f"Cannot connect to SQLite database {profile_path!r}: {err}")

    c = db.cursor()

    try:
        try:
            c.execute("SELECT strftime('%Y-%m-%d %H:%M:%f', t), label, t FROM mmpdb_times ORDER BY id")
        except sqlite3.DatabaseError as err:
            die(f"Cannot use SQLite database {profile_path!r}: {err}")
        except sqlite3.OperationalError as err:
            die(f"Cannot get times from {profile_path!r}: {err}")

        prev_t = None
        start_t = None
        click.echo(f"T-T0\tdelta_T\ttimestamp\tlabel")
        
        N_secs = 86400 # number of second in a day
        for datetime_str, label, t in c:
            if prev_t is None:
                prev_t = t
                start_t = t
            msg = "{:.3f}\t{:.3f}\t{}\t{}".format(
                (t-start_t)*N_secs, # elapsed since start
                (t-prev_t)*N_secs,  # delta since previous
                datetime_str,
                label,
                )
            click.echo(msg)
            prev_t = t
        
        c.close()
        
    finally:
        db.close()
        
