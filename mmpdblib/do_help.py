# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from __future__ import print_function, absolute_import

import sys
import argparse

def print_help_command(parser, args):
    sys.stdout.write(args.help_text.rstrip() + "\n\n")

#### mmpdb help-analysis

help_analysis_text = """\

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
  % mmpdb index test_data.fragments -o test_data.mmpa --properties test_data.csv
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

"""

    
help_analysis_text_old = """\

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
substructure in the variable part is a single fragment, while the
contant part contains one fragment for each cut.

The matched molecular pair indexing process finds all pairs which have
the same constant part, in order to define a transformation from one
variable part to another variable part. A 'rule' stores information
about a transformation, including a list of all the pairs for that
rule.

The 'rule environment' extendes the transformation to include
information about the local environment of the attachment points on
the constant part. The environment fingerprint is based on the RDKit
circular fingerprints for the attachment points. There is one rule
environment for each available radius. Larger radii correspond to more
specific environments. The 'rule environment statistics' table stores
information about the distribution of property changes for all of the
pairs which are contain the given rule and environment, with one table
for each property.


1) Fragment structures
======================

Use "smifrag" to see how a given SMILES is fragmented. Use "fragment"
to fragment all of the compounds in a SMILES file.

"mmpdb smifrag" is a diagnostic tool to help understand how a given
SMILES will be fragmented and to experiment with the different
fragmentation options. For example:

  % mmpdb smifrag c1ccccc1OC
                     |-------------  variable  -------------|       |---------------------  constant  --------------------
  #cuts | enum.label | #heavies | symm.class | smiles       | order | #heavies | symm.class | smiles           | with-H   
  ------+------------+----------+------------+--------------+-------+----------+------------+------------------+----------
    1   |     N      |    2     |      1     | [*]OC        |    0  |    6     |      1     | [*]c1ccccc1      | c1ccccc1 
    1   |     N      |    6     |      1     | [*]c1ccccc1  |    0  |    2     |      1     | [*]OC            | CO       
    2   |     N      |    1     |     11     | [*]O[*]      |   01  |    7     |     12     | [*]C.[*]c1ccccc1 | -        
    1   |     N      |    1     |      1     | [*]C         |    0  |    7     |      1     | [*]Oc1ccccc1     | Oc1ccccc1
    1   |     N      |    7     |      1     | [*]Oc1ccccc1 |    0  |    1     |      1     | [*]C             | C        


Use "mmpdb fragment" to fragment a SMILES file and produce a fragment
file for the MMP analysis. Start with the test data file named
'test_data.smi' containing the following structures:

Oc1ccccc1 phenol
Oc1ccccc1O catechol
Oc1ccccc1N 2-aminophenol
Oc1ccccc1Cl 2-chlorophenol
Nc1ccccc1N o-phenylenediamine
Nc1cc(O)ccc1N amidol
Oc1cc(O)ccc1O hydroxyquinol
Nc1ccccc1 phenylamine
C1CCCC1N cyclopentanol

  % mmpdb fragment test_data.smi -o test_data.fragments

Fragmentation can take a while. You can save time by asking the code
to reuse fragmentations from a previous run. If you do that then the
fragment command will reuse the old fragmentation parameters. (You
cannot override them with command-line options.). Here is an example:

  % mmpdb fragment data_file.smi -o new_data_file.fragments \\
         --cache old_data_file.fragments

The --cache option will greatly improve the fragment performance when
there are only a few changes from the previous run.

The fragmentation algorithm is configured to ignore structures which
are too big or have too many rotatable bonds. There are also options
which change where to make cuts and the number of cuts to make. Use
the '--help' option on each command for details.

Use "mmpdb help-smiles-format" for details about to parse different
variants of the SMILES file format.

2) Index the MMPA fragments to create a database
================================================

The "mmpa index" command indexes the output fragments from "mmpa
fragment" by their variable fragments, that is, it finds
fragmentations with the same R-groups and puts them together. Here's
an example:

  % mmpdb index test_data.fragments -o test_data.mmpdb

The output from this is a SQLite database.

If you have activity/property data and you do not want the database to
include structures where there is no data, then you can specify the
properties file as well:

  % mmpdb index test_data.fragments -o test_data.mmpa \\
       --properties test_data.csv

Use "mmpdb help-property-format" for property file format details.

For more help use "mmpdb index --help".

3) Add properties to a database
===============================

Use "mmpdb loadprops" to add or modify activity/property data in the
database. Here's an example property file named 'test_data.csv' with
molecular weight and melting point properties:

ID	MW	MP
phenol	94.1	41
catechol	110.1	105
2-aminophenol	109.1	174
2-chlorophenol	128.6	8
o-phenylenediamine	108.1	102
amidol	124.1	*
hydroxyquinol	126.1	140
phenylamine	93.1	-6
cyclopentanol	86.1	-19

The following loads the property data to the MMPDB database file
created in the previous section:

  % mmpdb loadprops -p test_data.csv test_data.mmpdb

Use "mmpdb help-property-format" for property file format details.

For more help use "mmpdb loadprops --help". Use "mmpdb list" to see
what properties are already loaded.

4) Identify possible transforms
===============================

Use "mmpdb transform" to transform an input structure using the rules
in a database. For each transformation it can estimate the effect on
any properties. The following looks at possible ways to transform
2-pyridone using the test dataset created in the previous section, and
predict the effect on the "MW" property (the output is reformatted for
clarity):

  % mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MW
  ID      SMILES MW_from_smiles MW_to_smiles  MW_radius  \\
   1  Clc1ccccn1         [*:1]O      [*:1]Cl          1
   2   Nc1ccccn1         [*:1]O       [*:1]N          1
   3    c1ccncc1         [*:1]O     [*:1][H]          1

                               MW_fingerprint  MW_rule_environment_id  \\
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     298
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     275
  tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     267

  MW_count  MW_avg  MW_std  MW_kurtosis  MW_skewness  MW_min  MW_q1  \\
         1    18.5     NaN          NaN          NaN    18.5   18.5
         3    -1.0     0.0          NaN          0.0    -1.0   -1.0
         4   -16.0     0.0          NaN          0.0   -16.0  -16.0

  MW_median  MW_q3  MW_max  MW_paired_t  MW_p_value
       18.5   18.5    18.5          NaN         NaN
       -1.0   -1.0    -1.0  100000000.0         NaN
      -16.0  -16.0   -16.0  100000000.0         NaN

This says that "c1cccnc1O" can be transformed to "Clc1ccccn1" using
the transformation [*:1]O>>[*:1]Cl (that is, replace the oxygen with a
chlorine). The best transformation match has a radius of 1, which
includes the aromatic carbon at the attachment point but not the
aromatic nitrogen which is one atom away.

There is only one pair for this transformation, and it predicts a shift
in molecular weight of 18.5. This makes sense as the [OH] is replaced
with a [Cl].

On the other hand, there are three pairs which transform it to
pyridine. The standard deviation of course 0 because it's a simple
molecular weight calculation. The 100000000.0 is the mmpdb way of
writing "positive infinity".

Melting point is more complicated. The following show that in the
transformation of 2-pyridone to pyridine there are still 3 matched
pairs and in this case the average shift is -93C with a standard
deviation of 76.727C:

  % mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MP
  ID      SMILES MP_from_smiles MP_to_smiles  MP_radius  \\
  1  Clc1ccccn1         [*:1]O      [*:1]Cl          1
  2   Nc1ccccn1         [*:1]O       [*:1]N          1
  3    c1ccncc1         [*:1]O     [*:1][H]          1

                               MP_fingerprint  MP_rule_environment_id  \\
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     298
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     275
 tLP3hvftAkp3EUY+MHSruGd0iZ/pu5nwnEwNA+NiAh8                     267

  MP_count  MP_avg  MP_std  MP_kurtosis  MP_skewness  MP_min   MP_q1  \\
        1 -97.000     NaN          NaN          NaN     -97  -97.00
        3 -16.667  75.235         -1.5     -0.33764     -72  -65.75
        3 -93.000  76.727         -1.5      0.32397    -180 -151.00

  MP_median  MP_q3  MP_max  MP_paired_t  MP_p_value
       -97 -97.00     -97          NaN         NaN
       -47  40.00      69       0.3837     0.73815
       -64 -42.25     -35       2.0994     0.17062

You might try enabling the --explain option to see why the algorithm
selected a given tranformation.

For more help use "mmpdb transform --help".

5) Use MMP to make a prediction
===============================

Use "mmpdb predict" to predict the property change in a transformation
from a given reference structure to a given query structure. Use this
when you want to limit the transform results when you know the
starting and ending structures. The following predicts the effect on
molecular weight in tranforming 2-pyridone to pyridone:

  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \\
            test_data.mmpdb --property MP
  predicted delta: -93 +/- 76.7268

This is the same MP_value and MP_std from the previous section using
'transform'.

  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \\
            test_data.mmpdb --property MP --value -41.6

I'll redo the calculation with the molecular weight property, and have
mmpdb do the trival calculation of adding the known weight to the
predicted delta:

  % mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \\
            test_data.mmpdb --property MW --value 95.1
  predicted delta: -16 predicted value: 79.1 +/- 0

You might try enabling the --explain option to see why the algorithm
selected a given transformation, or use --save-details to save list of
possible rules to the file 'pred_detail_rules.txt' and to save the
list of rule pairs to 'pred_detail_pairs.txt'.

For more help use "mmpdb predict --help".

"""
#### mmpdb help-admin

help_admin_text = """\

The administrative commands are:
  * list: describe what's inside of a database
  * loadprops: add or modify property information
  * smicat: show the structures in the database
  * propcat: show the properties in the database
  * drop_index: drop the database indices
  * create_index: (re)create the database indices

See the --help options for each command for more details.


"""


#### mmpdb help-smiles-format

help_smiles_format_text = """\

This help topic explains how the "--delimiter" and "--has-header"
options  of the "mmpa fragment" command affect SMILES parsing.

The mmpa code support the most common variants of a SMILES file. Every
SMILES file stores line-oriented records, with the SMILES in the first
field and the id (also called the title) in the second field. However,
there are differences in how to handle the first line of the file, and
in how to distinguish which is the second field. Some people use the
first line to store a header for each column in the file.

The classic Daylight SMILES file had no header line and interprets the
each line as a SMILES string followed by a whitespace followed by the
id/title. The id is the rest of the line, which means it may include
space and tabs. This is useful if you have identifiers with a space
in them, like IUPAC names or common names like "vitamin D".

A common variant is to treat the SMILES file as a CSV file, that is,
with at least two columns separated by a space, tab, or whitespace
character. Columns beyond the second column are ignored.

Use the "--delimiter" option to specify the delimiter type. The
available delimiter values are:

  "whitespace" (default) - CSV file with one or more whitespace
      characters as the delimiter
  "space" - CSV file with the space character as the delimiter
  "tab" - CSV file with the tab character as the delimiter
  "to-eol" - follow the Daylight rule where the id is
      everything past the first whitespace character
  "comma" - CSV file with a comma character as the delimiter
      (this is a very non-standard SMILES file format!)

The "native" delimiter is for chemfp compatibility. It is equivalent
to "whitespace", which matches the native RDKit parsing style.

Another variant takes the CSV format one step further, and lists
column headers as the first line of the file. The mmpa code will
ignore that line if you add the "--has-header" option.

Example command-line parameters:
   --delimiter whitespace       -- the default uses one or more whitespace
                                    as the column delimiter
   --delimiter to-eol           -- Daylight-style SMILES format
   --delimiter tab --has-header -- tab-delimited CSV, with a header line

"""


#### mmpdb help-fragments-format

# This appears to be out of date
help_fragments_format_text = """\

"*.fragments" format
====================

This is the output format for "mmpa fragment" and the input format for
"mmpa index". It is a line-oriented format, where each line is a JSON
record. This is known as a "JSON Lines" format; http://jsonlines.org/ .

Each line contains a record, represented as a JSON list. The first
term in the list is the record type. The five types are "VERSION",
"SOFTWARE", "OPTION", "RECORD", and "IGNORE".

The first line must contain a "VERSION" record, which specifies the
file format version. It must be the following:

  ["VERSION", "mmpdb-fragment/2"]

The second line is a "SOFTWARE" record, containing the name and
version of the program used to generate the fragments:

  ["SOFTWARE", "mmpdb-2.1"]

The "OPTION" records come next. These store the key/value pairs which
describe the fragmentation options. Each record is a 3-element
list. The second term is the key and the third is the value. The
following are the default mmpdb options:

  ["OPTION", "cut_smarts", "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]"]
  ["OPTION", "max_heavies", "100"]
  ["OPTION", "max_rotatable_bonds", "10"]
  ["OPTION", "method", "chiral"]
  ["OPTION", "num_cuts", "3"]
  ["OPTION", "rotatable_smarts", "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]"]
  ["OPTION", "salt_remover", "<default>"]

Afterwards are the RECORD and IGNORE records. If an input SMILES
record can be parsed and if it passes all of the fragment filters,
then the fragmentations are stored in a "RECORD" record, otherwise the
reason for the failure is stored in an "IGNORE" record.

Here is an annotated RECORD for 2-aminophenol ("Oc1ccccc1N"):

  ["RECORD",         #[0] "RECORD"
   "2-aminophenol",  #[1] the compound id from the SMILES file 
   "Oc1ccccc1N",     #[2] the input SMILES from the SMILES file
   8,                #[3] the number of heavies in the cleaned structure
   "Nc1ccccc1O",     #[4] the canonical SMILES of the cleaned structure
   [                 #[5] a list of fragmentations
    [1, "N", 1, "1", "[*]N", "0", 7, "1", "[*]c1ccccc1O", "Oc1ccccc1"],
    [1, "N", 1, "1", "[*]O", "0", 7, "1", "[*]c1ccccc1N", "Nc1ccccc1"],
    [1, "N", 7, "1", "[*]c1ccccc1N", "0", 1, "1", "[*]O", "O"],
    [1, "N", 7, "1", "[*]c1ccccc1O", "0", 1, "1", "[*]N", "N"],
    [2, "N", 6, "11", "[*]c1ccccc1[*]", "01", 2, "12", "[*]N.[*]O", null]]
   ]

Each fragmentation record is a 10 element list:
  [2,      #[0] number of cuts in the fragmentation
   "N",    #[1] chiral enumeration flag; "N" = No enumeration done
           #    "C" = constant up-enumeration, "V" = variable up-enumeration
   6,      #[2] number of heavy atoms in the variable part
   "11",   #[3] symmetry of the attachment points in the variable part
   "[*]c1ccccc1[*]",  #[4] the variable part as a canonical SMILES
   "01",   #[5] mapping from the variable parts attachment points to the constant
   2,      #[6] number of heavy atoms in the constant part
   "12",   #[7] symmetry of the attachment points in the constant part
   "[*]N.[*]O",  #[8] the constant part as a canonical SMILES
   null    #[9] for single cut structures, the canonical SMILES for the
           #     variable part attached to [H]; otherwise null.
  ]

The "IGNORE" records describe the input SMILES records which
were not fragmented. This can be because the SMILES could not
be parsed, or was too large, etc. The IGNORE record fields are:

  ["IGNORE",   #[0] "IGNORE"
   "WUNJUL",   #[1] the compound id from the SMIELS file
   "CCC(=O)CCCCCCCCC[C@@H]1CC[C@@H]([C@H](N1)CO)O",  #[2] the input SMILES
   "too many rotatable bonds"  #[3] the reason the structure was not fragmented
  ]

Some of the possible failure messages are: "invalid smiles",
"not enough heavy atoms", "too many heavy atoms", and
"too many rotatable bonds".

"""


help_properties_format_text = """\
This describes the file format used by the --properties option of
"mmpdb index" and "mmpdb loadprops" commands.

A property file contains information about the physical properties or
activity associated with a compound id. It is formatted as a data
table, where each line of the file is a row in the table. The first
line contains columns names. Each line contains exactly N fields, with
one field per column. The fields should be tab-separated. If the line
does not contain a tab character then it will be interpreted as
whitespace separated fields.

The first column contains compound identifiers. It must have the
column name "id", "ID", "Name", or "name". The remaining columns are
property columns, with the property name in the first row and property
values in the remaining rows. The property value for a given
identifier and property name must either be a number (something which
can be parsed as a floating point value; this includes integers) or
the symbol "*" to indicates that the value is missing.

Here is an example property file:

  ID MP CHR1 CHR2
  GEJYOJ 3 71 31.3
  ACIDUL 5 65 67.2
  KIXRIS 5 * *
  SOFWIV01 5 83 79.3

It defines three properties - MP, CHR1, and CHR2 - for four
identifiers. The MP value of GEJYOJ is 3, which is interpreted as
3.0. The CHR2 value of ACIDUL is interpted as 67.2. Compound KIXRIS
has the MP property 5.0 but the "*" indicates that its CHR1 and CHR2
properties are not known.

"""


def add_help_commands(subparsers):
    def add(subcommand, help, epilog):
        p = subparsers.add_parser(
            subcommand,
            help=help,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=epilog)
        p.set_defaults(command=print_help_command,
                       help_text=epilog,
                       subparser=p)
            
    add("help-analysis",
        "overview on how to use mmpdb for structure analysis",
        help_analysis_text)

    add("help-admin",
        "overview on how to use administor an mmpdb database",
        help_admin_text)

    add("help-smiles-format",
        "description of the SMILES parsing options",
        help_smiles_format_text)
        
    add("help-fragments-format",
        "description of the fragments file format",
        help_fragments_format_text)

    add("help-property-format",
        "description of the property file format",
        help_properties_format_text)
