CREATE INDEX fragmentation_on_record_id ON fragmentation(record_id);

-- Needed in fragdb_partition
CREATE INDEX fragmentation_on_constant_smiles ON fragmentation(constant_smiles);

CREATE UNIQUE INDEX record_on_title ON record(title);

CREATE UNIQUE INDEX error_record_on_title ON error_record(title);
