COPY index(constant_smiles,constant_symmetry_class,num_cuts,id,variable_symmetry_class,variable_smiles,attachment_order,enumeration_label)
FROM '/home/oriol/dev/mmpdb/new_db/"salida_single_cut.txt' DELIMITER '|' CSV HEADER;
