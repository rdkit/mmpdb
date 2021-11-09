CREATE INDEX fragmentation_on_record_id ON fragmentation(record_id);

CREATE INDEX record_on_title ON record(title);

CREATE INDEX error_record_on_title ON error_record(title);
