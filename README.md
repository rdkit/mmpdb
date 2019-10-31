# mmpdb 2.2-dev1 - matched molecular pair database generation and analysis


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

  "`mmpdb fragment`" -- fragment a SMILES file  
  "`mmpdb index`" -- find matched molecular pairs in a fragment file  

Use the "`--help`" option to get more information about any of the
commands. For example, "`mmpdb fragment --help`" will print the
command-line arguments, describe how they are used, and show
examples of use.

The subcommands starting with "help-" print additional information
about a given topic. The next few sections are the output from

```shell
   % mmpdb help-analysis
```

If you wish to experiment with a simple test set, use
tests/test_data.smi, with molecular weight and melting point
properties in tests/test_data.csv.


------------------


## Publication


An open-access publication describing this package has been 
published in the Journal of Chemical Information and Modeling:

A. Dalke, J. Hert, C. Kramer. mmpdb: An Open-Source Matched 
Molecular Pair Platform for Large Multiproperty Data Sets. *J. Chem. 
Inf. Model.*, **2018**, *58 (5)*, pp 902–910. 
https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173


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

Some terminology:
A fragmentation cuts 1, 2, or 3 non-ring bonds to
convert a structure into a "constant" part and a "variable" part. The
substructure in the variable part is a single fragment, and often
considered the R-groups, while the constant part contains one
fragment for each cut, and it often considered as containing the
core.

The matched molecular pair indexing process finds all pairs which have
the same constant part, in order to define a transformation from one
variable part to another variable part. A "rule" stores information
about a transformation, including a list of all the pairs for that
rule.

The "rule environment" extends the transformation to include
information about the local environment of the attachment points on
the constant part. The environment fingerprint is based on the RDKit
circular fingerprints for the attachment points. There is one rule
environment for each available radius. Larger radii correspond to more
specific environments. The "rule environment statistics" table stores
information about the distribution of property changes for all of the
pairs which contain the given rule and environment, with one table
for each property.



#### 1) Fragment structures


Use "`smifrag`" to see how a given SMILES is fragmented. Use "`fragment`"
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
file for the MMP analysis. Start with the test data file named
"test_data.smi" containing the following structures:

Oc1ccccc1 phenol  
Oc1ccccc1O catechol  
Oc1ccccc1N 2-aminophenol  
Oc1ccccc1Cl 2-chlorophenol  
Nc1ccccc1N o-phenylenediamine  
Nc1cc(O)ccc1N amidol  
Oc1cc(O)ccc1O hydroxyquinol  
Nc1ccccc1 phenylamine  
C1CCCC1N cyclopentanol  

```shell
  % mmpdb fragment test_data.smi -o test_data.fragments
```
Fragmentation can take a while. You can save time by asking the code
to reuse fragmentations from a previous run. If you do that then the
fragment command will reuse the old fragmentation parameters. (You
cannot override them with command-line options.). Here is an example:

```shell
  % mmpdb fragment data_file.smi -o new_data_file.fragments \ 
         --cache old_data_file.fragments
```

The "`--cache`" option will greatly improve the fragment performance when
there are only a few changes from the previous run.

The fragmentation algorithm is configured to ignore structures which
are too big or have too many rotatable bonds. There are also options
which change where to make cuts and the number of cuts to make. Use
the "`--help`" option on each command for details.

Use "`mmpdb help-smiles-format`" for details about to parse different
variants of the SMILES file format.

The "`--cut-smarts`" option sets the SMARTS pattern used to determine
which bonds to cut during fragmentation. Use "`--cut-rgroups`" or
"`--cut-rgroup-file`" to cut R-groups specified by fragment SMILES.

#### 2) Index the MMPA fragments to create a database


The "`mmpa index`" command indexes the output fragments from "`mmpa
fragment`" by their variable fragments, that is, it finds
fragmentations with the same R-groups and puts them together. Here's
an example:

```shell
  % mmpdb index test_data.fragments -o test_data.mmpdb
```
The output from this is a SQLite database.

If you have activity/property data and you do not want the database to
include structures where there is no data, then you can specify
the properties file as well:

```shell
  % mmpdb index test_data.fragments -o test_data.mmpdb --properties test_data.csv
```
Use "`mmpdb help-property-format`" for property file format details.

For more help use "`mmpdb index --help`".



#### 3) Add properties to a database


Use "`mmpdb loadprops`" to add or modify activity/property data in the
database. Here's an example property file named 'test_data.csv' with
molecular weight and melting point properties:

ID      MW      MP  
phenol  94.1    41  
catechol        110.1   105  
2-aminophenol   109.1   174  
2-chlorophenol  128.6   8  
o-phenylenediamine      108.1   102  
amidol  124.1   *  
hydroxyquinol   126.1   140  
phenylamine     93.1    -6  
cyclopentanol   86.1    -19  

The following loads the property data to the MMPDB database file
created in the previous section:

```shell
  % mmpdb loadprops -p test_data.csv test_data.mmpdb
```

Use "`mmpdb help-property-format`" for property file format details.

For more help use "`mmpdb loadprops --help`". Use "`mmpdb list`" to see
what properties are already loaded.



#### 4) Identify possible transforms


Use "`mmpdb transform`" to transform an input structure using the rules
in a database. For each transformation, it can estimate the effect on
any properties. The following looks at possible ways to transform
2-pyridone using the test dataset created in the previous section, and
predict the effect on the "MW" property (the output is reformatted for
clarity):

```shell
  % mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MW
  ID      SMILES MW_from_smiles MW_to_smiles  MW_radius  \ 
   1  Clc1ccccn1         [*:1]O      [*:1]Cl          1
   2   Nc1ccccn1         [*:1]O       [*:1]N          1
   3    c1ccncc1         [*:1]O     [*:1][H]          1

                               MW_fingerprint  MW_rule_environment_id  \ 
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     298
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     275
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     267

  MW_count  MW_avg  MW_std  MW_kurtosis  MW_skewness  MW_min  MW_q1  \ 
         1    18.5     NaN          NaN          NaN    18.5   18.5
         3    -1.0     0.0          NaN          0.0    -1.0   -1.0
         4   -16.0     0.0          NaN          0.0   -16.0  -16.0

  MW_median  MW_q3  MW_max  MW_paired_t  MW_p_value
       18.5   18.5    18.5          NaN         NaN
       -1.0   -1.0    -1.0  100000000.0         NaN
      -16.0  -16.0   -16.0  100000000.0         NaN
```

This says that "c1cccnc1O" can be transformed to "Clc1ccccn1" using
the transformation \[\*:1\]O>>\[\*:1\]Cl (that is, replace the oxygen with a
chlorine). The best transformation match has a radius of 1, which
includes the aromatic carbon at the attachment point but not the
aromatic nitrogen which is one atom away.

There is only one pair for this transformation, and it predicts a shift
in molecular weight of 18.5. This makes sense as the [OH] is replaced
with a [Cl].

On the other hand, there are three pairs which transform it to
pyridine. The standard deviation of course is 0 because it's a simple
molecular weight calculation. The 100000000.0 is the mmpdb way of
writing "positive infinity".

Melting point is more complicated. The following shows that in the
transformation of 2-pyridone to pyridine there are still 3 matched
pairs and in this case the average shift is -93C with a standard
deviation of 76.727C:

```shell
  % mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MP
  ID      SMILES MP_from_smiles MP_to_smiles  MP_radius  \ 
  1  Clc1ccccn1         [*:1]O      [*:1]Cl          1
  2   Nc1ccccn1         [*:1]O       [*:1]N          1
  3    c1ccncc1         [*:1]O     [*:1][H]          1

                               MP_fingerprint  MP_rule_environment_id  \ 
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     298
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     275
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     267

  MP_count  MP_avg  MP_std  MP_kurtosis  MP_skewness  MP_min   MP_q1  \ 
        1 -97.000     NaN          NaN          NaN     -97  -97.00
        3 -16.667  75.235         -1.5     -0.33764     -72  -65.75
        3 -93.000  76.727         -1.5      0.32397    -180 -151.00

  MP_median  MP_q3  MP_max  MP_paired_t  MP_p_value
       -97 -97.00     -97          NaN         NaN
       -47  40.00      69       0.3837     0.73815
       -64 -42.25     -35       2.0994     0.17062
```

You might try enabling the "`--explain`" option to see why the algorithm
selected a given tranformation.

For more help use "`mmpdb transform --help`".



#### 5) Use MMP to make a prediction


Use "`mmpdb predict`" to predict the property change in a transformation
from a given reference structure to a given query structure. Use this
when you want to limit the transform results when you know the
starting and ending structures. The following predicts the effect on
molecular weight in transforming 2-pyridone to pyridone:

```shell
  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \ 
            test_data.mmpdb --property MP
  predicted delta: -93 +/- 76.7268
```

This is the same MP_value and MP_std from the previous section using
'`transform`'.

```shell
  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \ 
            test_data.mmpdb --property MP --value -41.6
```

I'll redo the calculation with the molecular weight property, and have
mmpdb do the trival calculation of adding the known weight to the
predicted delta:

```shell
  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \ 
            test_data.mmpdb --property MW --value 95.1
  predicted delta: -16 predicted value: 79.1 +/- 0
```

You might try enabling the "`--explain`" option to see why the algorithm
selected a given transformation, or use "`--save-details`" to save the 
list of possible rules to the file 'pred_detail_rules.txt' and to save 
the list of rule pairs to "pred_detail_pairs.txt".


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


The mmpdb package is copyright 2015-2018 by F. Hoffmann-La
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

