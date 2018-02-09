select
	cs.smiles csmiles,
	c1.input_smiles pair1_smiles,
	c2.input_smiles pair2_smiles,
    rs1.smiles from_smiles,
	rs2.smiles to_smiles,
	re.radius,
	efp.fingerprint

from 
	compound c1, 
	compound c2, 
	pair p,
	constant_smiles cs,
	rule_environment re,
	rule r,
	rule_smiles rs1,
	rule_smiles rs2,
	environment_fingerprint efp
where 
	c1.id=p.compound1_id 
	and c2.id=p.compound2_id 
	and p.id=cs.id
	and p.rule_environment_id = re.id
	and re.rule_id=r.id
	and r.from_smiles_id = rs1.id
	and r.to_smiles_id = rs2.id
	and re.environment_fingerprint_id = efp.id