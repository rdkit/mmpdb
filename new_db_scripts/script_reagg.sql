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
from (select distinct * from index) i
group by constant_smiles,constant_symmetry_class,num_cuts;


CREATE INDEX ix_cuts
    ON public.index_agg USING btree
    (num_cuts ASC NULLS LAST)
    TABLESPACE pg_default;