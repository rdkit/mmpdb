-- mmpdb - matched molecular pair database generation and analysis
--
-- Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
--
-- Redistribution and use in source and binary forms, with or without
-- modification, are permitted provided that the following conditions are
-- met:
--
--    * Redistributions of source code must retain the above copyright
--      notice, this list of conditions and the following disclaimer.
--    * Redistributions in binary form must reproduce the above
--      copyright notice, this list of conditions and the following
--      disclaimer in the documentation and/or other materials provided
--      with the distribution.
--    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
--      its contributors may be used to endorse or promote products
--      derived from this software without specific prior written
--      permission.
--
-- THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
-- "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
-- LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
-- A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
-- HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
-- SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
-- LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
-- DATA, OR PROFITS, OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
-- THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
-- (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
-- OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--
-- END OF LICENSE

-- Create the tables for the mmpdb matched molecular pair database.

-- This schema is meant to support both SQLite and MySQL, but those
-- databases use slightly different SQL. This schema definition
-- contains two terms which must be replaced by text substution to
-- get the correct database-specific SQL. These are:

--   $PRIMARY_KEY$ - specifies the primary key type
--     SQLite: "INTEGER PRIMARY KEY"
--     MySQL:  "PRIMARY KEY"

--   $COLLATE$ - specifies a case-sensitive collation order
--     SQLite: ""
--     MySQL: "COLLATE latin1_bin"

-- In addition, to simplify processing, this file must only use
-- semicolons at the end of a SQL statement.

-- There's only one row in this table
CREATE TABLE dataset (
  id $PRIMARY_KEY$,
  mmpdb_version INTEGER NOT NULL,
  title VARCHAR(255) NOT NULL $COLLATE$,    -- human-visible label
  creation_date DATETIME NOT NULL,
  fragment_options VARCHAR(2000) NOT NULL $COLLATE$, -- the JSON-encoded options used to fragment the dataset
  index_options VARCHAR(2000) NOT NULL $COLLATE$, -- the JSON-encoded options used to index the dataset
  is_symmetric INTEGER NOT NULL,
  
  num_compounds INTEGER, -- used for "list"
  num_rules INTEGER,
  num_pairs INTEGER,
  num_rule_environments INTEGER,
  num_rule_environment_stats INTEGER
);



-- Each used input structure gets its own 'compound'
-- This table might not exist for rule-only data sets.

CREATE TABLE compound (
  id $PRIMARY_KEY$,
  public_id VARCHAR(255) NOT NULL $COLLATE$,    -- the public compound id (should this have a better name?)
  input_smiles VARCHAR(255) NOT NULL $COLLATE$, -- the input SMILES, before salt removal
  clean_smiles VARCHAR(255) NOT NULL $COLLATE$, -- the SMILES after salt removal
  clean_num_heavies INTEGER NOT NULL $COLLATE$  -- the number of heavies in the cleaned SMILES (needed?)
  );

-- Normalized property names (eg, "IC50" might be mapped to 3).

CREATE TABLE property_name (
  id $PRIMARY_KEY$,
  name VARCHAR(255) NOT NULL $COLLATE$
  );

-- Properties for each input compound

CREATE TABLE compound_property (
  id $PRIMARY_KEY$,
  compound_id INTEGER NOT NULL,
  property_name_id INTEGER NOT NULL,
  value REAL NOT NULL,
  FOREIGN KEY (compound_id) REFERENCES compound (id),
  FOREIGN KEY (property_name_id) REFERENCES property_name (id)
  );

-- -- A matched molecular pair rule

-- Normalized SMILES for the LHS or RHS of the rule transformation SMILES.
CREATE TABLE rule_smiles (
  id $PRIMARY_KEY$,
  smiles VARCHAR(255) NOT NULL $COLLATE$,
  num_heavies INTEGER
);



CREATE TABLE rule (
  id $PRIMARY_KEY$,

  -- The SMIRKS/transformation SMILES for this rule is:
  --   rule_smiles[id = from_smiles_id].smiles + ">>" +
  --   rule_smiles[id = to_smiles_id].smiles
  from_smiles_id INTEGER NOT NULL REFERENCES rule_smiles(id),
  to_smiles_id INTEGER NOT NULL REFERENCES rule_smiles(id)
  );

-- The "constant_part" (also called the context) is the substructure
-- which remains constant in the transformation. It typically has one
-- or more R-groups. It is omitted from the transformation A>>B.


-- A rule can have multiple rule environment, one per radius.

CREATE TABLE rule_environment (
 id INTEGER PRIMARY KEY,
 rule_id INTEGER REFERENCES rule(id),
 environment_fingerprint_id INTEGER REFERENCES environment_fingerprint(id),
 radius INTEGER
 -- should I keep track of the number of pairs?
 -- (the rule_env..._statistics "count" is the number of pairs with a given property)
 );

-- Table with normalized fingerprints for the rule_environment.

-- Fingerprints are based on the RDKit Morgan (ECFP-like) circular
-- fingerprints. They are only used for equality testing. To save
-- space, they are SHA2 hashed then base64 coded, so only 43 bytes
-- (344 bits) are needed.  They cannot be used for simliarity testing.

CREATE TABLE environment_fingerprint (
 id INTEGER PRIMARY KEY,
 fingerprint VARCHAR(43) NOT NULL $COLLATE$  -- the base64-encoded SHA2
 );


-- The pairs that belong to a rule_environment

CREATE TABLE pair (
  id $PRIMARY_KEY$,
  rule_environment_id INTEGER REFERENCES rule_environment (id) NOT NULL,
  compound1_id INTEGER NOT NULL REFERENCES compound (id),
  compound2_id INTEGER NOT NULL REFERENCES compound (id),
  constant_id INTEGER REFERENCES constant_smiles(id)
  );

CREATE TABLE constant_smiles (
  id $PRIMARY_KEY$,
  smiles VARCHAR(255)
  );

-- The aggregate property deltas for each rule environment

CREATE TABLE rule_environment_statistics (
  id $PRIMARY_KEY$,
  rule_environment_id INTEGER REFERENCES rule_environment (id),
  property_name_id INTEGER NOT NULL REFERENCES property_name (id),
  count INTEGER NOT NULL,
  avg REAL NOT NULL,
  std REAL,
  kurtosis REAL,
  skewness REAL,
  min REAL NOT NULL,
  q1 REAL NOT NULL,
  median REAL NOT NULL,
  q3 REAL NOT NULL,
  max REAL NOT NULL,
  paired_t REAL,
  p_value REAL
  );
