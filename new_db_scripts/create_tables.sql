drop table constant;

CREATE UNLOGGED TABLE public.constant (
	constant_smiles varchar(1000) NULL,
	constant_with_h_smiles varchar(1000) NULL
)
WITH (
	OIDS=FALSE
) ;
-- TODO PKey index -> UPSERT
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

drop table idrecord;
CREATE UNLOGGED TABLE public.idrecord (
	id varchar(1000) NOT NULL,
	input_smiles varchar(1000) NULL,
	num_normalized_heavies varchar(1000) NULL,
	normalized_smiles varchar(1000) NULL,
	PRIMARY KEY (id)
)
WITH (
	OIDS=FALSE
);
