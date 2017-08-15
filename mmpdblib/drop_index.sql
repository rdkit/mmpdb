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

DROP INDEX IF EXISTS rule_environment_rule_id;
DROP INDEX IF EXISTS rule_environment_environment_fingerprint_id;
DROP INDEX IF EXISTS rule_from_smiles_id;
DROP INDEX IF EXISTS rule_to_smiles_id;
DROP INDEX IF EXISTS rule_smiles_smiles;
DROP INDEX IF EXISTS environment_fingerprint_fingerprint;
DROP INDEX IF EXISTS rule_environment_statistics_rule_environment_and_property_name_ids;
DROP INDEX IF EXISTS rule_environment_statistics_count;

DROP INDEX IF EXISTS pair_rule_environment_id;
DROP INDEX IF EXISTS compound_property_compound_id_property_name_id;
