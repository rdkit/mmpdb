# CHANGELOG

## mmpdb 2.1 - (in development)

- RDKit 2018\_03\_1 changed the SMILES output so wildcard atoms are
  represented with a `*` instead of `[*]`. mmpdb works with SMILES
  strings at the syntax level. Parts of the code expected only the old
  RDKit behavior and crashed or gave the wrong output for the new
  behavior. For example, the fragmentation algorithm raised an
  AssertionError saying:

```
  File "mmpdblib/fragment_algorithm.py", line 368, in make_single_cut
    constant_smiles_with_H = replace_wildcard_with_H(constant_smiles)
  File "mmpdblib/fragment_algorithm.py", line 302, in replace_wildcard_with_H
    assert smiles.count("[*]") == 1, smiles
AssertionError: *O
```

  Version 2.1 now supports both ways of representing the wildcard atom
  and will work with RDKit releases from 2017 and 2018.

- The tests are now included as part of the distribution.

- mmpdb is available from PyPI using "pip install mmpdb"

- A preprint of 
  [our paper](https://chemrxiv.org/articles/mmpdb_An_Open_Source_Matched_Molecular_Pair_Platform_for_Large_Multi-Property_Datasets/5999375)
  is available from ChemRxiv. We received word today that
  [JCIM](https://pubs.acs.org/journal/jcisd8) has accepted it (after
  minor changes).

## mmpdb 2.0 - 2017-08-15

First public release of mmpdb, an open-source matched molecular pair
package designed to create and query MMP databases for big-pharma
sized ADMET datasets.
