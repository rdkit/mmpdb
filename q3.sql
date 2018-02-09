drop table mols_transformation;
create table mols_transformation as 
select distinct
	c1.input_smiles smiles_from,
	c2.input_smiles smiles_to,	        
         rs1.smiles smirks_from,
         rs2.smiles smirks_to,
          cs.smiles smiles_constant
from 
	compound c1, 
	compound c2, 
	pair p,	
	constant_smiles cs,
 	rule_environment re,
	rule r,
	rule_smiles rs1,
	rule_smiles rs2
where 
	c1.id=p.compound1_id 
	and c2.id=p.compound2_id
           and p.rule_environment_id = re.id
           and re.rule_id=r.id
           and cs.id = p.constant_id
           and r.from_smiles_id = rs1.id
           and r.to_smiles_id = rs2.id
order by  smiles_from,smiles_to

select smiles_from,count(*) from mols_transformation group by smiles_from order by count(*) desc

select count(*) from pair
select count(*) from compound