# CHANGELOG

## mmpdb 3.1 - 2023-11-28

Extended the "generate" command to handle 2-, and 3-cut transforms.

The `generate --explain` option now also explains why the search for a
matching query or variable part passes or fails. This proved useful in
determining that an expected fragmentation was instead being filtered.

The new `--min-heavies-total-const-frag N` fragmentation option
specifies the minimum number of heavy atoms allowed in the constant
part. The default value is 0.

Changed the fragdb schema version (in the "options" table) from 3 to 4
to support the new fragmentation option. Version 3 fragment databases
are still supported, by `min_heavies_total_const_frag` to 0.

Added two indices and a SQLite pragma for the page size. Roche reports
these improve analysis performance.

Fixed `--from` and `--to` support in proprulecat. These had been left
behind in the migration to click from argparse for command-line
processing.

## mmpdb 3.0 - 2023-5-31 

A large number of changes to merge three different development tracks
and add new features.

The "fragments" file format has been replaced with a SQLite-based
"fragdb" file format. This makes it much easier to develop tools to
work on fragment data sets instead of processing a JSON-Lines file.

New functionality to create an MMP data set in a distributed compute
environment. Some of the features are:

- split a SMILES file into a set of smaller SMILES files
- the default "fragment" file output is now based on the input name
- fragment files can be re-partitioned by constant fragments:
    - the "fragdb_constants" file generates fragment information
    - the "fragdb_partition" create re-partitioned fragdb files
- the default "index" file output is now based on the input name
- there are tools to merge fragdb and mmpdb files into one

As a result, mmpdb can now handle significantly larger data sets.

Added support for Postgres for direct index database creation. (The
new distributed compute tools require SQLite.)

Added a new "generate" command to apply 1-cut transforms to a
structure, using MMP rules as a playbook.

Replaced the SHA256-based Morgan fingerprint signature with a
canonical SMARTS representing the Morgan fingerprint environment. This
is difficult to understand or depict, so also include a "pseudo"
SMILES that can be parsed by RDKit (if sanitize is disabled) and
drawn. The new environment fingerprint also include the SMARTS of its
parent, that is, the SMARTS with a smaller radius.

Switched to 'click' for command-line parsing, removed the vendered
version of the peewee ORM, and switched to a modern "pyproject.toml"
project configuration with a setup.cfg which declares its dependencies.


## mmpdb 2.2-dev (the GitHub development track)

The `fragment` and `smifrag` commands now support options for
supervised fragmentation based on a specified set of R-group SMILES to
use for the fragmentation. Multiple SMILES can be specified on the
command-line using the `--cut-rgroup` argument, with one SMILES per
argument, or using the `--cut-rgroup-file` argument with the name of
an R-group file. The file must be formatted with one R-group SMILES
per line.

All SMILES strings must contain a single wildcard atom ("*"), which
indicates the attachment point. The wildcard atom must contain only
one single bond to the rest of the R-group, and cannot contain charge,
hydrogens, isotope, or other properties.

The SMILES strings are converted into a SMARTS pattern which matchs
the SMILES exactly (each atom must have the same valence and hydrogen
count, and the bond types must match). These SMARTS patterns are then
merged into a single recusive SMARTS with two terms: the wildcard
atom, single bonded to a recursive SMARTS term for each of the SMILES
strings.

The new `rgroup2smarts` command can be used to process the R-group
SMILES into SMARTS, as a way to examine the conversion process and
verify that it works.

## mmpdb 2.2 - 2019-01-11

  This minor release contains improvements that help reducing the 
  database size. Many transformations and associated statistics inside
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
