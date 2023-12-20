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

-- NOTE: There is configuration information in three files!
-- 1) fragment_types.py -- the data types
-- 2) fragment_schema.sql -- (this file) defines the SQL schema
-- 3) fragment_db.py -- defines the mapping from SQL to the data types


CREATE TABLE options (
	id INTEGER NOT NULL,
	version INTEGER,
	cut_smarts VARCHAR(1000),
	max_heavies INTEGER,
	max_rotatable_bonds INTEGER,
	method VARCHAR(20),
	num_cuts INTEGER,
	rotatable_smarts VARCHAR(1000),
	salt_remover VARCHAR(200),
	min_heavies_per_const_frag INTEGER,
	min_heavies_total_const_frag INTEGER,
        max_up_enumerations INTEGER,
	PRIMARY KEY (id)
);

CREATE TABLE error_record (
	id INTEGER NOT NULL,
	title VARCHAR(100) NOT NULL,
	input_smiles VARCHAR(300) NOT NULL,
	errmsg VARCHAR(100),
	PRIMARY KEY (id)
);

-- Unfortunately, this 'record' does not use the same colum names as
-- 'compound' in the mmpdb schema.

CREATE TABLE record (
	id INTEGER NOT NULL,
	title VARCHAR(50) NOT NULL,
	input_smiles VARCHAR(400) NOT NULL,
	num_normalized_heavies INTEGER,
	normalized_smiles VARCHAR(350) NOT NULL,
	PRIMARY KEY (id)
);

CREATE TABLE fragmentation (
	id INTEGER NOT NULL,
	record_id INTEGER,
	num_cuts INTEGER,
	enumeration_label VARCHAR(1) NOT NULL,
	variable_num_heavies INTEGER,
	variable_symmetry_class VARCHAR(3) NOT NULL,
	variable_smiles VARCHAR(350) NOT NULL,
	attachment_order VARCHAR(3) NOT NULL,
	constant_num_heavies INTEGER,
	constant_symmetry_class VARCHAR(3) NOT NULL,
	constant_smiles VARCHAR(350) NOT NULL,
	constant_with_H_smiles VARCHAR(350),
	PRIMARY KEY (id),
	FOREIGN KEY(record_id) REFERENCES record (id)
);
