-- mmpdb - matched molecular pair database generation and analysis
--
-- Copyright (c) 2021, Andrew Dalke Scientific, AB
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


-- Version 3.0 switched to a SQLite database to store the fragments.
-- Earlier versions used JSON-Lines.
-- The SQLite database improves I/O time, reduces memory use, and
--     simplifies the development of fragment analysis tools.
CREATE TABLE config (
  version INTEGER,
  cut_smarts TEXT,
  max_heavies INTEGER,
  max_rotatable_bonds INTEGER,
  method STRING,
  num_cuts INTEGER,
  rotatable_smarts STRING,
  salt_remover STRING
);

CREATE TABLE error_record (
  id INTEGER PRIMARY KEY,
  title STRING,
  input_smiles STRING,
  errmsg STRING
);

CREATE TABLE record (
  id INTEGER PRIMARY KEY,
  title TEXT,
  input_smiles TEXT NOT NULL,
  num_normalized_heavies INTEGER,
  normalized_smiles TEXT NOT NULL
  );

CREATE TABLE fragmentation (
  id INTEGER PRIMARY KEY,
  record_id INTEGER REFERENCES record(id),
  num_cuts INTEGER,
  enumeration_label TEXT NOT NULL,
  variable_num_heavies INTEGER,
  variable_symmetry_class TEXT NOT NULL,
  variable_smiles TEXT NOT NULL,
  attachment_order TEXT NOT NULL,
  constant_num_heavies INTEGER,
  constant_symmetry_class INTEGER,
  constant_smiles TEXT NOT NULL,
  constant_with_H_smiles TEXT
);

