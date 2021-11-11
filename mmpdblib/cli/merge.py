import os
import click

from .click_utils import (
    command,
    die,
    add_multiple_databases_parameters,
    )

MERGE_CREATE_INDEX_SQL = """
CREATE UNIQUE INDEX idx_compound ON compound(public_id);
CREATE UNIQUE INDEX idx_rule_smiles ON rule_smiles(smiles);
CREATE UNIQUE INDEX idx_rule ON rule (from_smiles_id, to_smiles_id);
CREATE UNIQUE INDEX idx_env_fp ON environment_fingerprint (smarts, parent_smarts);
CREATE UNIQUE INDEX idx_rule_env ON rule_environment (rule_id, environment_fingerprint_id, radius);
CREATE UNIQUE INDEX idx_constant_smiles ON constant_smiles(smiles);

CREATE TEMPORARY TABLE rule_map (
  old_id INTEGER PRIMARY KEY,
  new_id INTEGER
);

CREATE TEMPORARY TABLE rule_env_map (
  old_id INTEGER PRIMARY KEY,
  new_id INTEGER
);

"""

MERGE_DROP_INDEX_SQL = """
DROP INDEX idx_compound;
DROP INDEX idx_rule_smiles;
DROP INDEX idx_rule;
DROP INDEX idx_env_fp;
DROP INDEX idx_rule_env;
DROP INDEX idx_constant_smiles;
"""

MERGE_DATABASE_SQL = """
-- This expects the database to import to be attached as 'old'
-- using something like:
--   attach database "subset.000.mmpdb" as old

-- Step 1: Copy over the compound table

-- Must have UNIQUE INDEX ON compound (public_id)
-- The merge code ensured duplicate entries have the same clean_smiles.

INSERT OR IGNORE INTO compound (public_id, input_smiles, clean_smiles, clean_num_heavies)
 SELECT public_id, input_smiles, clean_smiles, clean_num_heavies
   FROM old.compound
        ;

-- Step 2: Copy over the rule_smiles

-- Must have UNIQUE INDEX ON rule_smiles (smiles)

INSERT OR IGNORE INTO rule_smiles (smiles, num_heavies)
 SELECT smiles, num_heavies
   FROM old.rule_smiles
        ;

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

-- Step 4. Table mapping old rule id to new rule id

-- This simplifies steps 6 and 7

INSERT INTO rule_map (old_id, new_id) 
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


-- Step 5: Merge environment_fingerprint

-- Must have UNIQUE INDEX ON environment_fingerprint (smarts, parent_smarts)
-- (I haven't checked if I only need to use the smarts.)

INSERT OR IGNORE INTO environment_fingerprint (smarts, pseudosmiles, parent_smarts)
  SELECT smarts, pseudosmiles, parent_smarts
    FROM old.environment_fingerprint
         ;

-- Step 6: Merge rule_environment

-- Must have UNIQUE INDEX ON rule_environment (rule_id, environment_fingerprint_id, radius)


INSERT OR IGNORE INTO rule_environment (rule_id, environment_fingerprint_id, radius)
  SELECT rule_map.new_id, new_env_fp.id, old_rule_env.radius
    FROM old.rule_environment AS old_rule_env,
         rule_map,
         old.environment_fingerprint AS old_env_fp,
         environment_fingerprint AS new_env_fp
   WHERE old_rule_env.rule_id = rule_map.old_id
     AND old_rule_env.environment_fingerprint_id = old_env_fp.id
     AND (    old_env_fp.smarts = new_env_fp.smarts
          AND (   old_env_fp.parent_smarts = new_env_fp.parent_smarts
               OR (old_env_fp.parent_smarts IS NULL AND new_env_fp.parent_smarts IS NULL)))
         ;

-- Step 7: Table mapping old rule environment to new rule environment


INSERT INTO rule_env_map (old_id, new_id)
  SELECT old_rule_env.id, new_rule_env.id
    FROM old.rule_environment AS old_rule_env,
         rule_environment AS new_rule_env,
         rule_map,
         old.environment_fingerprint AS old_env_fp,
         environment_fingerprint AS new_env_fp
   WHERE old_rule_env.rule_id = rule_map.old_id
     AND new_rule_env.rule_id = rule_map.new_id
     AND old_rule_env.environment_fingerprint_id = old_env_fp.id
     AND (    old_env_fp.smarts = new_env_fp.smarts
          AND (   old_env_fp.parent_smarts = new_env_fp.parent_smarts 
               OR (old_env_fp.parent_smarts IS NULL AND new_env_fp.parent_smarts IS NULL)))
     AND (    new_rule_env.environment_fingerprint_id = new_env_fp.id
          AND new_rule_env.radius = old_rule_env.radius
          AND new_rule_env.rule_id = rule_map.new_id)
         ;


-- Step 8: Merge constant_smiles

-- requires a unique index

INSERT OR IGNORE INTO constant_smiles (smiles)
   SELECT smiles
     FROM old.constant_smiles
          ;

-- Step 9: Merge pairs

-- These are all unique because the constant_smiles wasn't seen earlier

INSERT INTO pair (rule_environment_id, compound1_id, compound2_id, constant_id) 
 SELECT rule_env_map.new_id, new_compound1.id, new_compound2.id, new_constant.id
   FROM old.pair AS old_pair,
        rule_env_map,
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

DELETE FROM rule_env_map;
DELETE FROM rule_map;

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
        ):
    import sqlite3
    from .. import index_writers, reporters
    
    try:
        writer = index_writers.open_sqlite_index_writer(output_filename, title)
    except IOError as err:
        die(f"Cannot open SQLite file {output_filename!r}: {err}")
    except sqlite3.OperationalError as err:
        die(f"Unexpected problem initializing mmpdb file {output_filename!r}: {err}")

    writer.start(fragment_options, index_options)

    # Create the index
    writer.end(reporters.Quiet())
    writer.conn.execute("COMMIT")

    try:
        return sqlite3.connect(output_filename)
    except sqlite3.OperationalError as err:
        die(f"Unexpected problem re-opening mmpdb file {output_filename!r}: {err}")

        
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

@add_multiple_databases_parameters()
@click.pass_obj
def merge(
        reporter,
        databases_options,
        title,
        output_filename,
        ):
    """merge multiple mmpdb databases

    Each DATABASE must be an mmpdb file.

    Properties are ignored because merging them is not meaningful.
    """
    import sqlite3
    from .. import dbutils
    from .. import index_writers
    from .. import schema
    
    
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
        for database in databases:
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
                        )
                    output_c = output_db.cursor()
                    try:
                        schema._execute_sql(output_c, MERGE_CREATE_INDEX_SQL)
                    except Exception:
                        output_c.close()
                        output_db.close()
                        output_c = output_db = None
                        raise

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
            try:
                output_db.execute("ATTACH DATABASE ? AS old", (database,))
            except sqlite3.OperationalError as err:
                die(f"Cannot attach {database!r} to {output_file!r}: {err}")

            try:
                # Check for duplicate constants
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

                # Merge
                schema._execute_sql(output_c, MERGE_DATABASE_SQL)

            except:
                raise
            else:
                output_c.execute("DETACH DATABASE old")
                

    finally:
        if output_c is not None:
            schema._execute_sql(output_c, MERGE_DROP_INDEX_SQL)
            index_writers.update_counts(output_c)
            output_c.execute("ANALYZE")
            output_c.execute("COMMIT")
            output_c.close()
            output_db.close()

