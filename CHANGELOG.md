# CHANGELOG

## mmpdb 2.2 - 2019-11-01

  This minor release contains improvements that help reducing the 
  database size. MAny transformations and associated statistics inside
  the database are unlikely to ever be used, since there are other 
  transformations that will yield the same compounds. To accomplish this, 
  three new options have been introduced:

- max-radius: The maximum radius can now be set on the command line 
  during indexing. This was hardcoded in previous versions.

- smallest-transformation-only: Some transformation scan be reduced to 
  smaller transformations, for example p-Fluoro-phenyl >> p-Chloro-phenyl
  to Fluoro >> Chloro. If this flag is set during indexing, reducible 
  transformations will not be written to the database. Note that this only
  setting reduces the number of transformations for a given pair. It does
  not completely remove a pair.

- min-heavies-per-const-frag: For double- and triple-cuts, fragmentations
  are created where the constant parts can be very small down to a single
  atom. For example, if the fragmentation algorithm can generate a single-cut 
  where the variable fragment is p-Fluorophenyl, it will also generate a
  double cut with phenyl as variable fragment and F as one of the constant
  pieces. The 'min-heavies-per-const-frag' option can be used during 
  fragmentation to eliminate multiple cuts where one of the constant fragments
  is very small. If this is set to 3 or 4, double- and triple cuts only happen
  at positions that would be considered real scaffolds changes in the middle 
  of molecules. Note that in principle this option only reduces the number 
  of pairs for a given transformation, effectively removing multiple-cuts 
  where possible. There may be edge cases where pairs are completely removed,
  because the single cut transfers too many atoms. If you use this option, you 
  may want to adjust the --max-variable-heavies option during indexing.

  The last '--smallest-transformation-only' and '--min-heavies-per-const 3' 
  options together typically reduce the database size by ~ 70%.


## mmpdb 2.1 - 2018-04-27

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
