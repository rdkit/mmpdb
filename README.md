# mmpdb 2.0 - matched molecular pair database generation and analysis


## Synopsis


A package to identify matched molecular pairs and use them to predict
property changes.


------------------


## Requirements


The package has been tested on both Python 2.7 and Python 3.6.

You will need a copy of the RDKit cheminformatics toolkit, available
from http://rdkit.org/ . Apart from other standard scientific python 
libraries like scipy and numpy, this is the only required third-party
dependency for normal operation, though several optional third-party
packages may be used if available.

  - The matched molecular pairs are stored in a SQLite database. The
APSW module from https://github.com/rogerbinns/apsw gives slightly
better analysis performance than Python's built-in SQLite module.

  - The fragment file is in JSON Lines format (see http://jsonlines.org/ ).
The ujson (https://github.com/esnme/ultrajson) and slightly slower
cjson (https://github.com/AGProjects/python-cjson ) are both about
25% faster than Python 2.7's built-in 'json' module.

  - The "`--memory`" option in the index command requires the psutil
module (see https://pypi.python.org/pypi/psutil/5.2.2 ) to get memory
use information.


------------------


## How to run the program and get help


The package includes a command-line program named "mmpdb". This
support many subcommands. For examples:

  "mmpdb fragment" -- fragment a SMILES file
  "mmpdb index" -- find matched molecular pairs in a fragment file

Use the "`--help`" option to get more information about any of the
commands. For example, "mmpdb fragment --help" will print the
command-line arguments, describe how they are used, and show
examples of use.

The subcommands starting with 'help-' print additional information
about a given topic. The next few sections are the output from

```shell
   % mmpdb help-analysis
```

If you wish to experiment with a simple test set, use
tests/test_data.smi, with molecular weight and melting point
properties in tests/test_data.csv.


------------------


## Background


The overall process is:

  1) Fragment structures in a SMILES file, to produce fragments.

  2) Index the fragments to produces matched molecular pairs.
     (you might include property information at this point)

  3) Load property information.

  4) Find transforms for a given structure; and/or

  5) Predict a property for a structure given the known
     property for another structure

Some terminology. A fragmentation cuts 1, 2, or 3 non-ring bonds to
convert a structure into a "constant" part and a "variable" part. The
substructure in the variable part is a single fragment, and often
considered the core structure, while the constant part contains one
fragment for each cut, and it often considered as containing the
R-groups.

The matched molecular pair indexing process finds all pairs which have
the same constant part, in order to define a transformation from one
variable part to another variable part. A 'rule' stores information
about a transformation, including a list of all the pairs for that
rule.

The 'rule environment' extends the transformation to include
information about the local environment of the attachment points on
the constant part. The environment fingerprint is based on the RDKit
circular fingerprints for the attachment points. There is one rule
environment for each available radius. Larger radii correspond to more
specific environments. The 'rule environment statistics' table stores
information about the distribution of property changes for all of the
pairs which are contain the given rule and environment, with one table
for each property.


#### 1) Fragment structures


Use "smifrag" to see how a given SMILES is fragmented. Use "fragment"
to fragment all of the compounds in a SMILES file.

"`mmpdb smifrag`" is a diagnostic tool to help understand how a given
SMILES will be fragmented and to experiment with the different
fragmentation options. For example:

```shell
  % mmpdb smifrag 'c1ccccc1OC'
                     |-------------  variable  -------------|       |---------------------  constant  --------------------
  #cuts | enum.label | #heavies | symm.class | smiles       | order | #heavies | symm.class | smiles           | with-H   
  ------+------------+----------+------------+--------------+-------+----------+------------+------------------+----------
    1   |     N      |    2     |      1     | [*]OC        |    0  |    6     |      1     | [*]c1ccccc1      | c1ccccc1 
    1   |     N      |    6     |      1     | [*]c1ccccc1  |    0  |    2     |      1     | [*]OC            | CO       
    2   |     N      |    1     |     11     | [*]O[*]      |   01  |    7     |     12     | [*]C.[*]c1ccccc1 | -        
    1   |     N      |    1     |      1     | [*]C         |    0  |    7     |      1     | [*]Oc1ccccc1     | Oc1ccccc1
    1   |     N      |    7     |      1     | [*]Oc1ccccc1 |    0  |    1     |      1     | [*]C             | C        
```

Use "`mmpdb fragment`" to fragment a SMILES file and produce a fragment
file for the MMP analysis. For example:

```shell
  % mmpdb fragment csd.smi -o csd.fragments
```
Fragmentation can take a while. You can save time by asking the code
to reuse fragmentations from a previous run. If you do that then the
fragment command will reuse the old fragmentation parameters. (You
cannot override them with command-line options.). Here is an example:

```shell
  % mmpdb fragment datafile.smi -o new_datafile.fragments \
         --cache old_datafile.fragments
```

The "`--cache`" option will greatly improve the fragment performance when
there are only a few changes from the previous run.

The fragmentation algorithm is configured to ignore structures which
are too big or have too many rotatable bonds. There are also options
which change where to make cuts and the number of cuts to make. Use
the "`--help`" option on each command for details.

Use "mmpdb help-smiles-format" for details about to parse different
variants of the SMILES file format.


#### 2) Index the MMPA fragments to create a database


The "mmpa index" command indexes the output fragments from "mmpa
fragment" by their variable fragments, that is, it finds
fragmentations with the same R-groups and puts them together. Here's
an example:

```shell
  % mmpdb index csd.fragments -o csd.mmpdb
```
The output from this is a SQLite database.

If you have activity/property data and you do not the database to
include structures where there is no data, then you can specify
the properties file as well:

```shell
  % mmpdb index csd.fragments -o csd.mmpa --properties MP.csv
```
Use "`mmpdb help-property-format`" for property file format details.

For more help use "`mmpdb index --help`".


#### 3) Add properties to a database


Use "`mmpdb loadprops`" to add or modify activity/property data in the
database. The following adds melting point data in CSV format.

```shell
  % mmpdb loadprops -p MP.csv csd.mmpdb
```

Use "`mmpdb help-property-format`" for property file format details.

For more help use "`mmpdb loadprops --help`". Use "`mmpdb list`" to see
what properties are already loaded.


#### 4) Identify possible transforms


Use "`mmpdb transform`" to transform an input structure using the rules
in a database. For each transformation, it can estimate the effect on
any properties. The following looks at possible ways to transform
diphenyl ether where there are at least 100 known matched pairs based
on that transform.

```shell
  % mmpdb transform --smiles 'c1ccccc1Oc1ccccc1' csd.mmpdb --min-pairs 100

   ID      SMILES MP_from_smiles MP_to_smiles  MP_radius
    1  COc1ccccc1  [*:1]c1ccccc1       [*:1]C          0
    2   Oc1ccccc1  [*:1]c1ccccc1     [*:1][H]          0

                                MP_fingerprint  MP_rule_environment_id
   59SlQURkWt98BOD1VlKTGRkiqFDbG6JVkeTJ3ex3bOA                      90
   59SlQURkWt98BOD1VlKTGRkiqFDbG6JVkeTJ3ex3bOA                     956

   MP_count  MP_avg  MP_std  MP_kurtosis  MP_skewness  MP_min  MP_q1
        106 -23.443  36.987      1.56380      0.65077    -159    -44
        112 -11.696  39.857      0.18439      0.14109    -135    -38

   MP_median  MP_q3  MP_max  MP_paired_t    MP_p_value
       -20.0   -3.0      49       6.5256  2.447100e-09
       -11.5   13.5      99       3.1057  2.410100e-03
```

For more help use "`mmpdb transform --help`".


#### 5) Use MMP to make a prediction


Use "`mmpdb predict`" to predict the property change in a transformation
from a given reference structure to a given query structure. The
following estimates that the melting point of diphenyl ether will
increase by 5 degrees, with a standard deviation of 18 degrees, if the
oxygen is replaced by a sulfur.

```shell
  % mmpdb predict --reference 'c1ccccc1Oc1ccccc1' --smiles 'c1ccccc1Sc1ccccc1' csd.mmpdb -p MP
  predicted delta: +4.91667 +/- 18.0477
```

For more help use "mmpdb predict --help".


------------------


## History and Acknowledgements


The project started as a fork of the matched molecular pair program
'mmpa' written by Jameed Hussain, then at GlaxoSmithKline Research &
Development Ltd.. Many thanks to them for contributing the code to the
RDKit project under a free software license.

Since then it has gone through two rewrites. Major changes to the
first version included:
  - performance improvements,

  - support for property prediction

  - environmental fingerprints
  
That version supported both MySQL and SQLite, and used the third-party
"peewee.py" and "playhouse" code to help with for database
portability. Many thanks to Charlies Leifer for that software.

The second version dropped MySQL support but added APSW support, which
was already available in the peewee/playhouse modules. The major goals
in version 2 were:

  - better support for chiral structures

  - canonical variable fragments, so the transforms are canonical
    on both the left-hand and right-hand sides. (Previously only
    the entire transform was canonical.)


------------------


## Copyright


The mmpdb package is copyright 2015-2017 by F. Hoffmann-La
Roche Ltd and distributed under the 3-clause BSD license. See [LICENSE](LICENSE)
for details.


------------------


## License information


The software derives from software which is copyright 2012-2013 by
GlaxoSmithKline Research & Development Ltd., and distributed under the
3-clause BSD license. To the best of our knowledge, mmpdb does not contain any
of the mmpa original source code. We thank the authors for releasing this
package and include their license in the credits. See [LICENSE](LICENSE) for details.

The file fileio.py originates from [chemfp](http://chemfp.com) and is therefore
copyright by Andrew Dalke Scientific AB under the MIT license. See
[LICENSE](LICENSE) for details. Modifications to this file are covered under
the mmpdb license.

The files peewee.py and playhouse/\*.py are copyright 2010 by Charles
Leifer and distributed under the MIT license. See [LICENSE](LICENSE) for details.

