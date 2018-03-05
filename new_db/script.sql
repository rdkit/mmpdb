drop table constant;

CREATE UNLOGGED TABLE public.constant (
	constant_smiles varchar(1000) NULL,
	constant_with_h_smiles varchar(1000) NULL
)
WITH (
	OIDS=FALSE
) ;

drop table index;
CREATE UNLOGGED TABLE public.index (
	constant_smiles varchar(1000) NULL,
	constant_symmetry_class varchar(1000) NULL,
	num_cuts varchar NULL,
	id varchar(1000) NULL,
	variable_symmetry_class varchar(1000) NULL,
	variable_smiles varchar(1000) NULL,
	attachment_order varchar(1000) NULL,
	enumeration_label varchar(1000) NULL
)
WITH (
	OIDS=FALSE
) ;

drop table main;
CREATE UNLOGGED TABLE public.main (
	normalized_smiles varchar(1000) NULL,
	id varchar(1000) NULL
)
WITH (
	OIDS=FALSE
) ;


truncate table index;
COPY index(constant_smiles,constant_symmetry_class,num_cuts,id,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label)
FROM '/home/oriol/dev/mmpdb/new_db/salida.txt' DELIMITER '|' CSV HEADER;

truncate table constant;
COPY constant(constant_smiles,constant_with_H_smiles)
FROM '/home/oriol/dev/mmpdb/new_db/salida_cte.txt' DELIMITER '|' CSV HEADER;

truncate table main;
COPY main(normalized_smiles,id)
FROM '/home/oriol/dev/mmpdb/new_db/salida_global.txt' DELIMITER '|' CSV HEADER;

CREATE INDEX ix_index
    ON public.index USING btree
    (constant_smiles ASC NULLS LAST, constant_symmetry_class ASC NULLS LAST, num_cuts ASC NULLS LAST)
    TABLESPACE pg_default;

CREATE INDEX ix_ct
    ON public.constant USING btree
    (constant_smiles ASC NULLS LAST)
    TABLESPACE pg_default;

create UNLOGGED table constant_unique
as select distinct constant_smiles,constant_with_h_smiles from constant;

CREATE INDEX ix_ct2
    ON public.constant_unique USING btree
    (constant_smiles ASC NULLS LAST)
    TABLESPACE pg_default;

CREATE INDEX ix_main
    ON public.main USING btree
    (normalized_smiles ASC NULLS LAST)
    TABLESPACE pg_default;

drop type index_value_type;
CREATE TYPE index_value_type AS (
  id character varying(1000) ,
  variable_symmetry_class character varying(1000) ,
    variable_smiles character varying(1000),
    attachment_order character varying(1000),
    enumeration_label character varying(1000)
);

drop table index_agg;
create UNLOGGED table index_agg as
SELECT constant_smiles,constant_symmetry_class,num_cuts,
  jsonb_agg(
      cast(ROW(id ,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label) as index_value_type)
  ) as variable_part
from index
group by constant_smiles,constant_symmetry_class,num_cuts;

CREATE INDEX ix_cuts
    ON public.index_agg USING btree
    (num_cuts ASC NULLS LAST)
    TABLESPACE pg_default;