# mmpdb 3.1 - matched molecular pair database generation and analysis


## Synopsis


A package to identify matched molecular pairs and use them to predict
property changes and generate new molecular structures.


------------------

## Installation

mmpdb 3.1 must be installed before use. (Earlier versions of mmpdb
could be run in-place, in the top-level directory.) This will also
ensure that the SciPy, peewee, and click packages are installed.

To install from PyPI using
[pip](https://pip.pypa.io/en/stable/user_guide/), which comes with
Python:

On macOS and other Unix-like systems:
```
python -m pip install mmpdb
```
On Windows:
```
py -m pip install mmpdb
```

If you are using a virtual environment (eg, with
[venv](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment)
or
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
then the installer will place a copy of the relevant files into the
virtual environment's Python package directory, and place the
`mmpdb` command-line driver in your path.

If you are not using a virtual environment (and you should be using a
virtual environment) then the installer will place the relevant files
into Python's system directory. If you do not have write permissions
to that directory, or only want to install it for your personal use,
then add the `--user` flag at the end of the install command.

To install from source directory, go to the top-level directory then
do:

On macOS and other Unix-like systems:
```
python -m pip install .
```

On Windows:
```
py -m pip install .
```

If you plan to modify the source code and want the installation to use
the local directory rather than make a copy into the package
directory, then you need an "editable" installation,

On macOS and other Unix-like systems:
```
python -m pip install -e .
```

On Windows:
```
py -m pip install -e .
```

(Assuming the dependencies are installed then it is possible to use
mmpdb without installation by going to the top-level directory and
using `python -m mmpdblib` or, for Windows, `py -m mmpdlib`. This is
not recommended.)

## Requirements

The package has been tested on Python 3.9 and 3.10. It should work
under newer versions of Python.

You will need a copy of the RDKit cheminformatics toolkit, available
from http://rdkit.org/ , which in turn requires NumPy. You will also
need SciPy, peewee, and click. The latter three are listed as
dependencies in setup.cfg and should be installed automatically.

Optional components you may find useful are:

  - The matched molecular pairs may instead by be stored in a Postgres
database. These were tested using the
[psycopg2](https://www.psycopg.org/) adapter. See `mmpdb
help-postgres` for more information.

 - The "`--memory`" option in the index command requires the
[psutil](https://pypi.python.org/pypi/psutil/) module to get memory
use information.

NOTE: mmpdb 2 used a JSON-Lines format for the fragment files, and
suggested an optional package with faster JSON parsing. mmpdb 3 no
longer uses this format.


------------------


## How to run the program and get help


The package includes a command-line program named "mmpdb". This
support many subcommands. For examples:

* "`mmpdb fragment`" -- fragment a SMILES file

* "`mmpdb index`" -- find matched molecular pairs in a fragment file

Use the "`--help`" option to get more information about any of the
commands. For example, "`mmpdb fragment --help`" will print the
command-line arguments, describe how they are used, and show
examples of use.

The subcommands starting with "help-" print additional information
about a given topic. Much of the text of this README come from the
output of

```shell
 % mmpdb help-analysis 
 % mmpdb help-distributed 
```

If you wish to experiment with a simple test set, use
`tests/test_data.smi`, with molecular weight and melting point
properties in `tests/test_data.csv`.


------------------


## Publication


An open-access publication describing this package has been 
published in the Journal of Chemical Information and Modeling:

A. Dalke, J. Hert, C. Kramer. mmpdb: An Open-Source Matched 
Molecular Pair Platform for Large Multiproperty Data Sets. *J. Chem. 
Inf. Model.*, **2018**, *58 (5)*, pp 902–910. 
https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173

For more about the methods to scale mmpdb to larger datasets and
generate new molecules, see:

M. Awale, J. Hert, L. Guasch, S. Riniker, C. Kramer.
The Playbooks of Medicinal Chemistry Design Moves. *J. Chem. 
Inf. Model.*, **2021**,  *61 (2)*, pp 729-742.
https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c01143


------------------


## Background

The overall process is:

1) Fragment structures in a SMILES file, to produce fragments.

2) Index the fragments to produces matched molecular pairs. (you might include
property information at this point)

3) Load property information.

4) Find transforms for a given structure; and/or

5) Predict a property for a structure given the known    property for another
structure; and/or

6) Apply 1-cut rules to generate new structures from a given    structure.

Some terminology:

A fragmentation cuts 1, 2, or 3 non-ring bonds to convert a structure into a
"constant" part and a "variable" part. The substructure in the variable part
is a single fragment, and often considered the R-groups, while the constant
part contains one fragment for each cut, and it often considered as containing
the core.

The matched molecular pair indexing process finds all pairs which have the
same constant part, in order to define a transformation from one variable part
to another variable part. A "rule" stores information about a transformation,
including a list of all the pairs for that rule.

The "rule environment" extends the transformation to include information about
the local environment of the attachment points on the constant part. The
environment fingerprint is based on the RDKit circular fingerprints for the
attachment points, expressed as a canonical SMARTS pattern, and alternatively,
as a "pseudo"-SMILES string, which is a bit less precise but easier to
understand and visualize.

The fingerprint SMARTS pattern describes the Morgan circular fingerprint
invariants around the attachment points. Here's a 2-cut example split across
three lines:

```
[#0;X1;H0;+0;!R:1]-[#6;X4;H1;+0;R](-[#6;X4;H2;+0;R])-[#6;X4;H2;+0;R].
[#0;X1;H0;+0;!R:2]-[#7;X3;H0;+0;R](-[#6;X4;H2;+0;R])-[#6;X4;H2;+0;R].
[#0;X1;H0;+0;!R:3]-[#6;X3;H0;+0;R](:[#6;X3;H1;+0;R]):[#6;X3;H1;+0;R]
```

The SMARTS modifiers, like "H0" to require no hydrogens, are needed to match
the Morgan invariants but are quite the eye-full. The psuedosmiles alternative
is:

```
[*:1]-[CH](-[CH2](~*))-[CH2](~*).
[*:2]-[N](-[CH2](~*))-[CH2](~*).
[*:3]-[c](:[cH](~*)):[cH](~*)
```

This can be processed by RDKit, if sanitization is disabled, and turned into
an image.

CAUTION! The "`(~*)`" terms are used to represent the SMARTS connectivity
terms "X<digit>", but they do not necessarily all represent distinct atoms!

There is one rule environment for each available radius. Larger radii
correspond to more specific environments. The "rule environment statistics"
table stores information about the distribution of property changes for all of
the pairs which contain the given rule and environment, with one table for
each property.

### 1) Fragment structures

Use "`smifrag`" to see how a given SMILES is fragmented. Use "`fragment`" to
fragment all of the compounds in a SMILES file.

"`mmpdb smifrag`" is a diagnostic tool to help understand how a given SMILES
will be fragmented and to experiment with the different fragmentation options.
For example:

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

Use "`mmpdb fragment`" to fragment a SMILES file and produce a fragment file
for the MMP analysis. Start with the test data file named "test_data.smi"
containing the following structures:

```text
Oc1ccccc1 phenol  
Oc1ccccc1O catechol  
Oc1ccccc1N 2-aminophenol  
Oc1ccccc1Cl 2-chlorophenol  
Nc1ccccc1N o-phenylenediamine  
Nc1cc(O)ccc1N amidol  
Oc1cc(O)ccc1O hydroxyquinol  
Nc1ccccc1 phenylamine  
C1CCCC1N cyclopentanol  
```

then run the following command generate a fragment database.

```shell
% mmpdb fragment test_data.smi -o test_data.fragdb
```

Fragmentation can take a while. You can save time by asking the code to reuse
fragmentations from a previous run. If you do that then the fragment command
will reuse the old fragmentation parameters. (You cannot override them with
command-line options.). Here is an example:

```shell
% mmpdb fragment data_file.smi -o new_data_file.fragdb \
       --cache old_data_file.fragdb
```

The "`--cache`" option will greatly improve the fragment performance when
there are only a few changes from the previous run.

The fragmentation algorithm is configured to ignore structures which are too
big or have too many rotatable bonds. There are also options which change
where to make cuts and the number of cuts to make. Use the "`--help`" option
on each command for details.

Use "`mmpdb help-smiles-format`" for details about to parse different variants
of the SMILES file format.

### 2) Index the MMPA fragments to create a database

The "`mmpa index`" command indexes the output fragments from "`mmpa fragment`"
by their variable fragments, that is, it finds fragmentations with the same
R-groups and puts them together. Here's an example:

```shell
% mmpdb index test_data.fragdb -o test_data.mmpdb
```

The output from this is a SQLite database.

If you have activity/property data and you do not want the database to include
structures where there is no data, then you can specify the properties file as
well:

```shell
% mmpdb index test_data.fragdb -o test_data.mmpdb --properties test_data.csv
```

Use "`mmpdb help-property-format`" for more details about the property file
format.

For more help use "`mmpdb index --help`".

### 3) Add properties to a database

Use "`mmpdb loadprops`" to add or modify activity/property data in the
database. Here's an example property file named 'test_data.csv' with molecular
weight and melting point properties:

```text
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
```

The following loads the property data to the MMPDB database file created in
the previous section:

```shell
% mmpdb loadprops -p test_data.csv test_data.mmpdb
Using dataset: MMPs from 'test_data.fragdb'
Reading properties from 'tests/test_data.csv'
Read 2 properties for 9 compounds from 'tests/test_data.csv'
Imported 9 'MW' records (9 new, 0 updated).
Imported 8 'MP' records (8 new, 0 updated).
Number of rule statistics added: 533 updated: 0 deleted: 0
Loaded all properties and re-computed all rule statistics.
```

Use "`mmpdb help-property-format`" for more details about the property file
format.

For more help use "`mmpdb loadprops --help`". Use "`mmpdb list`" to see what
properties are already loaded.

### 4) Identify possible transforms

Use "`mmpdb transform`" to transform an input structure using the rules in a
database. For each transformation, it can estimate the effect on any
properties. The following looks at possible ways to transform 2-pyridone using
the test dataset created in the previous section, and predict the effect on
the "MW" property (the output is reformatted for clarity):

```shell
% mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MW
ID     SMILES    MW_from_smiles    MW_to_smiles    MW_radius
1    Clc1ccccn1     [*:1]O          [*:1]Cl           1
2     Nc1ccccn1     [*:1]O          [*:1]N            1
3     c1ccncc1      [*:1]O          [*:1][H]          1

      MW_smarts                        MW_pseudosmiles    MW_rule_environment_id 
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]    [*:1]-[#6](~*)(~*)    299
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]    [*:1]-[#6](~*)(~*)    276
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]    [*:1]-[#6](~*)(~*)    268

MW_count    MW_avg    MW_std    MW_kurtosis    MW_skewness
    1        18.5
    3        -1         0            0
    4       -16         0            0

MW_min  MW_q1  MW_median  MW_q3  MW_max  MW_paired_t    MW_p_value
 18.5    18.5   18.5      18.5    18.5
 -1      -1     -1        -1      -1      1e+08    
-16     -16    -16       -16     -16      1e+08 
```

This says that "c1cccnc1O" can be transformed to "Clc1ccccn1" using the
transformation \[\*:1\]O>>\[\*:1\]Cl (that is, replace the oxygen with a
chlorine). The best transformation match has a radius of 1, which includes the
aromatic carbon at the attachment point but not the aromatic nitrogen which is
one atom away.

There is only one pair for this transformation, and it predicts a shift in
molecular weight of 18.5. This makes sense as the [OH] is replaced with a
[Cl].

On the other hand, there are three pairs which transform it to pyridine. The
standard deviation of course is 0 because it's a simple molecular weight
calculation. The 1e+08.0 is the mmpdb way of writing "positive infinity".

Melting point is more complicated. The following shows that in the
transformation of 2-pyridone to pyridine there are still 3 matched pairs and
in this case the average shift is -93C with a standard deviation of 76.727C:

```shell
% mmpdb transform --smiles 'c1cccnc1O' test_data.mmpdb --property MP
ID   SMILES    MP_from_smiles   MP_to_smiles   MP_radius   
1   Clc1ccccn1    [*:1]O           [*:1]Cl        1
2    Nc1ccccn1    [*:1]O           [*:1]N         1
3    c1ccncc1     [*:1]O           [*:1][H]       1

MP_smarts                            MP_pseudosmiles     MP_rule_environment_id
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]   [*:1]-[#6](~*)(~*)   299
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]   [*:1]-[#6](~*)(~*)   276
[#0;X1;H0;+0;!R:1]-[#6;X3;H0;+0;R]   [*:1]-[#6](~*)(~*)   268

MP_count   MP_avg   MP_std   MP_kurtosis   MP_skewness   
   1       -97            
   3       -16.667   75.235     -1.5        -0.33764   
   3       -93       76.727     -1.5        -0.32397   

MP_min   MP_q1   MP_median   MP_q3   MP_max   MP_paired_t   MP_p_value
 -97      -97       -97       -97     -97      
 -72      -65.75    -47        40      69        0.3837      0.73815
-180     -151       -64       -42.25  -35       -2.0994      0.17062
```

You might try enabling the "`--explain`" option to see why the algorithm
selected a given tranformation.

For more help use "`mmpdb transform --help`".

### 5) Use MMP to make a prediction

Use "`mmpdb predict`" to predict the property change in a transformation from
a given reference structure to a given query structure. Use this when you want
to limit the transform results when you know the starting and ending
structures. The following predicts the effect on molecular weight in
transforming 2-pyridone to pyridone:

```shell
% mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \
          test_data.mmpdb --property MP
predicted delta: -93 +/- 76.7268
```

This is the same MP_value and MP_std from the previous section using
'`transform`'.

The reference value may also be included in the calulation, to give a
predicted value.

```shell
% mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \
          test_data.mmpdb --property MP --value -41.6
predicted delta: -93 predicted value: -134.6 +/- 76.7268
```

I'll redo the calculation with the molecular weight property, and have mmpdb
do the trival calculation of adding the known weight to the predicted delta:

```shell
% mmpdb predict --smiles 'c1cccnc1' --reference 'c1cccnc1O' \
          test_data.mmpdb --property MW --value 95.1
predicted delta: -16 predicted value: 79.1 +/- 0
```

You might try enabling the "`--explain`" option to see why the algorithm
selected a given transformation, or use "`--save-details`" to save the  list
of possible rules to the file `pred_detail_rules.txt` and to save  the list of
rule pairs to `pred_detail_pairs.txt`.

### 6) Use MMP to generate new structures

The rules in a MMP database give a sort of "playbook" about the
transformations which might be explored in medicinal chemistry. These rule can
be applied to a given structure to generate new related structures, following
a method related to the transform command but ignoring any property
information. Here's an example using the default radius of 0, which means the
environment fingerprint is ignored. (The columns have been re-formatted for
the documentation.)

```shell
% mmpdb generate --smiles 'c1ccccc1C(O)C' test_data.mmpdb
start             constant  from_smiles  to_smiles          r  pseudosmiles  final
CC(O)c1ccccc1  *C(C)c1ccccc1  [*:1]O    [*:1][H]            0  [*:1](~*)  CCc1ccccc1
CC(O)c1ccccc1  *C(C)c1ccccc1  [*:1]O    [*:1]N              0  [*:1](~*)  CC(N)c1ccccc1
CC(O)c1ccccc1  *C(C)c1ccccc1  [*:1]O    [*:1]Cl             0  [*:1](~*)  CC(Cl)c1ccccc1
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccccc1O      0  [*:1](~*)  CC(O)c1ccccc1O
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccccc1N      0  [*:1](~*)  CC(O)c1ccccc1N
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1cc(O)ccc1N   0  [*:1](~*)  CC(O)c1cc(O)ccc1N
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccc(O)cc1N   0  [*:1](~*)  CC(O)c1ccc(O)cc1N
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]C1CCCC1        0  [*:1](~*)  CC(O)C1CCCC1
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccccc1Cl     0  [*:1](~*)  CC(O)c1ccccc1Cl
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccc(N)c(N)c1 0  [*:1](~*)  CC(O)c1ccc(N)c(N)c1
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1cc(O)ccc1O   0  [*:1](~*)  CC(O)c1cc(O)ccc1O
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccc(O)c(O)c1 0  [*:1](~*)  CC(O)c1ccc(O)c(O)c1
CC(O)c1ccccc1  *C(C)O     [*:1]c1ccccc1 [*:1]c1ccc(O)cc1O   0  [*:1](~*)  CC(O)c1ccc(O)cc1O

#pairs  pair_from_id  pair_from_smiles  pair_to_id  pair_to_smiles
4       2-aminophenol  Nc1ccccc1O     phenylamine        Nc1ccccc1
3       phenol         Oc1ccccc1      phenylamine        Nc1ccccc1
1       catechol       Oc1ccccc1O     2-chlorophenol     Oc1ccccc1Cl
2       phenylamine    Nc1ccccc1      2-aminophenol      Nc1ccccc1O
2       phenylamine    Nc1ccccc1      o-phenylenediamine Nc1ccccc1N
1       phenylamine    Nc1ccccc1      amidol             Nc1ccc(O)cc1N
1       phenylamine    Nc1ccccc1      amidol             Nc1ccc(O)cc1N
1       phenylamine    Nc1ccccc1      cyclopentanol      NC1CCCC1
1       phenol         Oc1ccccc1      2-chlorophenol     Oc1ccccc1Cl
1       phenol         Oc1ccccc1      amidol             Nc1ccc(O)cc1N
1       phenol         Oc1ccccc1      hydroxyquinol      Oc1ccc(O)c(O)c1
1       phenol         Oc1ccccc1      hydroxyquinol      Oc1ccc(O)c(O)c1
1       phenol         Oc1ccccc1      hydroxyquinol      Oc1ccc(O)c(O)c1
```

The second half the output shows the number of known pairs for the given rule
environment (use `--min-pairs N` to require at least N pairs), and gives a
representative pair from the dataset.

In the above example, all of the fragmentations in the specified `--smiles`
are used. Alternatively, you may specify `--smiles` and one of `--constant` or
`--query` to use that specific fragmentation, or use `--constant` and
`--query` (without `--smiles`) to specify the exact pair.

There is also an option to generate `--subqueries`. This generates all of the
unique 1-cut fragmentations of the query, and uses them as additional queries.
I'll use the `--constant` to specify the phynol group, leaving the
aminomethanol available as the query. I'll use `--subqueries` to include
fragments of the query. I'll limit the output `--columns` to the start and
final SMILES structures, and the number of pairs. I'll use `--explain` to
display debug information, and finally, I'll use `--no-header` to make the
output a bit less complicated:

```shell
% mmpdb generate --smiles 'c1ccccc1C(O)N' --constant '*c1ccccc1' test_data.mmpdb \
     --subqueries --columns start,final,#pairs --explain --no-header
Number of subqueries: 4
Subqueries are: ['*CN', '*CO', '*N', '*O']
Using constant SMILES *c1ccccc1 with radius 0.
Environment SMARTS: [#0;X1;H0;+0;!R:1] pseudoSMILES: [*:1](~*)
Number of matching environment rules: 42
Query SMILES [*:1]C(N)O is not a rule_smiles in the database.
Query SMILES [*:1]CN is not a rule_smiles in the database.
Query SMILES [*:1]CO is not a rule_smiles in the database.
Nc1ccccc1     Oc1ccccc1       3
Nc1ccccc1     c1ccccc1        2
Nc1ccccc1     Clc1ccccc1      1
Number of rules for [*:1]N: 3
Oc1ccccc1     c1ccccc1        4
Oc1ccccc1     Nc1ccccc1       3
Oc1ccccc1     Clc1ccccc1      1
Number of rules for [*:1]O: 3
```

## Distributed computing

These commands enable MMP generation on a distributed compute cluster, rather
than a single machine.

NOTE: This method does not support properties, and you must use the
SQLite- based "mmpdb" files, not Postgres databases. The
[Postgres wiki](https://wiki.postgresql.org/wiki/Converting_from_other_Databases_to_PostgreSQL)
mentions [pgloader](https://github.com/dimitri/pgloader) as a possible
tool to have Postgres load a SQLite database.

These examples assume you work in a queueing environment with a shared
file system, and a queueing system which lets you submit a command and
a list of filenames, to enqueue the command once for each filename.

This documentation will use the command 'qsub' as a wrapper around [GNU
Parallel](https://www.gnu.org/software/parallel/):

```shell
alias qsub="parallel --no-notice -j 1 --max-procs 4"
```

This alias suppresses the request to cite GNU parallel in scientific papers,
and has it process one filename at a time, with at most 4 processes in
parallel.

I'll pass the filenames to process via stdin, like this example:

```shell
% ls /etc/passwd ~/.bashrc | qsub wc
       2       5      88 /Users/dalke/.bashrc
     120     322    7630 /etc/passwd
```

This output shows that `wc` received only a single filename because with two
filenames it also shows a 'total' line.

```shell
% wc /etc/passwd ~/.bashrc
     120     322    7630 /etc/passwd
       2       5      88 /Users/dalke/.bashrc
     122     327    7718 total
```

### Distributed fragmentation generation

NOTE: This method can also be used to process larger data sets on a single
machine because the `mmpdb merge` step uses less memory than the `mmpdb
index`.

The `fragment` command supports multi-processing with the `-j` flag, which
scales to about 4 or 8 processors. For larger data sets you can break the
SMILES dataset into multiple files, fragment each file indepenently, then
merge the results.

These steps are:

* smi_split - split the SMILES file into smaller files
* fragment - fragment the each smaller SMILES file into its own fragb file.
* fragdb_merge - merge the smaller fragdb files together.

#### Use smi_split to create N smaller SMILES files

I'll start with a SMILES file containing a header and 20267 SMILES lines:

```shell
% head -3 ChEMBL_CYP3A4_hERG.smi
SMILES  CMPD_CHEMBLID
[2H]C([2H])([2H])Oc1cc(ncc1C#N)C(O)CN2CCN(C[C@H](O)c3ccc4C(=O)OCc4c3C)CC2       CHEMBL3612928
[2H]C([2H])(N[C@H]1C[S+]([O-])C[C@@H](Cc2cc(F)c(N)c(O[C@H](COC)C(F)(F)F)c2)[C@@H]1O)c3cccc(c3)C(C)(C)C  CHEMBL2425617
% wc -l ChEMBL_CYP3A4_hERG.smi
   20268 ChEMBL_CYP3A4_hERG.smi
```

By default the "smi_split" command splits a SMILES file into 10 files. (Use
`-n` or `--num-files` to change the number of files, or use `--num-records` to
have N records per file.)

```shell
% mmpdb smi_split ChEMBL_CYP3A4_hERG.smi
Created 10 SMILES files containing 20268 SMILES records.
```

That "20268 SMILES record" shows that all 20268 lines were used to generate
SMILES records, which is a mistake as it includes the header line. I'll re-do
the command with `--has-header` to have it skip the header:

```shell
% mmpdb smi_split ChEMBL_CYP3A4_hERG.smi --has-header
Created 10 SMILES files containing 20267 SMILES records.
```

By default this generates files which look like:

```shell
% ls -l ChEMBL_CYP3A4_hERG.*.smi
-rw-r--r--  1 dalke  admin  141307 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0000.smi
-rw-r--r--  1 dalke  admin  152002 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0001.smi
-rw-r--r--  1 dalke  admin  127397 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0002.smi
-rw-r--r--  1 dalke  admin  137930 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0003.smi
-rw-r--r--  1 dalke  admin  130585 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0004.smi
-rw-r--r--  1 dalke  admin  150072 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0005.smi
-rw-r--r--  1 dalke  admin  139620 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0006.smi
-rw-r--r--  1 dalke  admin  133347 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0007.smi
-rw-r--r--  1 dalke  admin  131310 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0008.smi
-rw-r--r--  1 dalke  admin  129344 Feb 10 15:10 ChEMBL_CYP3A4_hERG.0009.smi
```

The output filenames are determined by the `--template` option, which defaults
to `{prefix}.{i:04}.smi`, where `i` is the output file index. See `smi_split
--help` for details.

#### Fragment the SMILES files

These files can be fragmented in parallel:

```shell
% ls ChEMBL_CYP3A4_hERG.*.smi | qsub mmpdb fragment -j 1
```

I used the `-j 1` flag to have `mmpdb fragment` use only a single thread,
otherwise each of the four fragment commands will use 4 threads even though my
laptop only has 4 cores. You should adjust the value to match the resources
available on your compute node.

The `parallel` command doesn't forward output until the program is done, so it
takes a while to see messages like:

```
Using 'ChEMBL_CYP3A4_hERG.0002.fragdb' as the default --output file.
Fragmented record 249/2026 (12.3%)[15:04:16] Conflicting single bond
directions around double bond at index 5.
[15:04:16]   BondStereo set to STEREONONE and single bond directions set to NONE.
```

If no `-o`/`--output` is specified, the `fragment` command uses a named based
on the input name, for example, if the input file is
`ChEMBL_CYP3A4_hERG.0002.smi` then the default output file is
`ChEMBL_CYP3A4_hERG.0002.mmpdb`.

#### Merge the fragment files

NOTE: This step is only needed if you want to use the merged file as a
`--cache` for new fragmentation. The `fragdb_constants` and `fragdb_partition`
commands can work directly on the un-merged fragdb files.

About 28 minutes later I have 10 fragdb files:

```shell
% ls -l ChEMBL_CYP3A4_hERG.*.fragdb
-rw-r--r--  1 dalke  admin  17862656 Feb 10 15:17 ChEMBL_CYP3A4_hERG.0000.fragdb
-rw-r--r--  1 dalke  admin  38285312 Feb 10 15:27 ChEMBL_CYP3A4_hERG.0001.fragdb
-rw-r--r--  1 dalke  admin  15024128 Feb 10 15:16 ChEMBL_CYP3A4_hERG.0002.fragdb
-rw-r--r--  1 dalke  admin  15929344 Feb 10 15:16 ChEMBL_CYP3A4_hERG.0003.fragdb
-rw-r--r--  1 dalke  admin  18063360 Feb 10 15:23 ChEMBL_CYP3A4_hERG.0004.fragdb
-rw-r--r--  1 dalke  admin  20586496 Feb 10 15:24 ChEMBL_CYP3A4_hERG.0005.fragdb
-rw-r--r--  1 dalke  admin  24911872 Feb 10 15:26 ChEMBL_CYP3A4_hERG.0006.fragdb
-rw-r--r--  1 dalke  admin  16875520 Feb 10 15:28 ChEMBL_CYP3A4_hERG.0007.fragdb
-rw-r--r--  1 dalke  admin  12451840 Feb 10 15:28 ChEMBL_CYP3A4_hERG.0008.fragdb
-rw-r--r--  1 dalke  admin  11010048 Feb 10 15:29 ChEMBL_CYP3A4_hERG.0009.fragdb
```

I'll merge these with the `fragdb_merge` command:

```shell
% mmpdb fragdb_merge ChEMBL_CYP3A4_hERG.*.fragdb -o ChEMBL_CYP3A4_hERG.fragdb
Merge complete. #files: 10 #records: 18759 #error records: 1501
```

This took about 4 seconds.

#### Use the merged fragment file as cache

The merged file can be used a a cache file for future fragmentations, such as:

```shell
% ls ChEMBL_CYP3A4_hERG.*.smi | \
    qsub mmpdb fragment --cache ChEMBL_CYP3A4_hERG.fragdb -j 1
```

This re-build using cache takes about 20 seconds.

# Distributed indexing

The `mmpdb index` command is single-threaded. It's possible to parallelize
indexing by partitioning the fragments with the same constant SMILES into
their own fragdb data sets, indexing those files, then merging the results
back into a full MMP database.

Note: the merge command can only be used to merge MMP databases with distinct
constants. It cannot be used to merge arbitrary MMP databases.

Note: the MMP database only stores aggregate information about pair
properties, and the aggregate values cannot be meaningfully merged, so the
merge command will ignore any properties in the database.

#### Partitioning on all constants

The `mmpdb fragdb_partition` command splits one or more fragment databases
into N smaller files. All of the fragmentations with the same constant are in
the same file.

NOTE: the fragdb files from the `fragment` command have a slightly different
structure than the ones from the `partition` command. The fragment fragdb
files only contain the input records that were fragmented. Each partition
fragdb file contains *all* of the input records from the input fragment
file(s). This is needed to handle 1-cut hydrogen matched molecular pairs.

If you specify multiple fragdb files then by default the results are put into
files matching the template "partition.{i:04d}.fragdb", as in the following:

```shell
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb
Analyzed 'ChEMBL_CYP3A4_hERG.0000.fragdb': #constants: 48895 #fragmentations: 109087
Analyzed 'ChEMBL_CYP3A4_hERG.0001.fragdb': #constants: 70915 #fragmentations: 212777
Analyzed 'ChEMBL_CYP3A4_hERG.0002.fragdb': #constants: 52370 #fragmentations: 100594
Analyzed 'ChEMBL_CYP3A4_hERG.0003.fragdb': #constants: 49021 #fragmentations: 103350
Analyzed 'ChEMBL_CYP3A4_hERG.0004.fragdb': #constants: 52318 #fragmentations: 112930
Analyzed 'ChEMBL_CYP3A4_hERG.0005.fragdb': #constants: 55977 #fragmentations: 123463
Analyzed 'ChEMBL_CYP3A4_hERG.0006.fragdb': #constants: 64083 #fragmentations: 164259
Analyzed 'ChEMBL_CYP3A4_hERG.0007.fragdb': #constants: 51605 #fragmentations: 114113
Analyzed 'ChEMBL_CYP3A4_hERG.0008.fragdb': #constants: 44149 #fragmentations: 80613
Analyzed 'ChEMBL_CYP3A4_hERG.0009.fragdb': #constants: 35889 #fragmentations: 69029
Analyzed 10 databases. Found #constants: 467865 #fragmentations: 1190215
Exporting 1 constants to 'partition.0000.fragdb' (#1/10, weight: 334589647)
Exporting 1 constants to 'partition.0001.fragdb' (#2/10, weight: 270409141)
Exporting 1 constants to 'partition.0002.fragdb' (#3/10, weight: 225664391)
Exporting 1 constants to 'partition.0003.fragdb' (#4/10, weight: 117895691)
Exporting 77977 constants to 'partition.0004.fragdb' (#5/10, weight: 52836587)
Exporting 77978 constants to 'partition.0005.fragdb' (#6/10, weight: 52836587)
Exporting 77975 constants to 'partition.0006.fragdb' (#7/10, weight: 52836586)
Exporting 77976 constants to 'partition.0007.fragdb' (#8/10, weight: 52836586)
Exporting 77977 constants to 'partition.0008.fragdb' (#9/10, weight: 52836586)
Exporting 77978 constants to 'partition.0009.fragdb' (#10/10, weight: 52836586)
```

The command's `--template` option lets you specify how to generate the output
filenames.

Why are there so few constants in first files and so many in the other? And
what are the "weight"s?

I'll use the `fragdb_constants` command to show the distinct constants in each
file and the number of occurrences.

```shell
% mmpdb fragdb_constants partition.0000.fragdb
constant        N
*C      25869
```

That's a lot of methyls (25,869 to be precise).

The indexing command does `N*(N-1)/2` indexing comparisions, plus a 1-cut
hydrogen match, so the cost estimate for the methyls is `25869*(25869-1)/2+1 =
334589647`, which is the `weight` value listed above.

I'll next list the three most common and least constants in
ChEMBL_CYP3A4_hERG.0004.fragdb:

```shell
% mmpdb fragdb_constants partition.0004.fragdb --limit 3
constant        N
*C.*C.*OC       7076
*C.*Cl  4388
*C.*C.*CC       3261
% mmpdb fragdb_constants partition.0004.fragdb | tail -3
*n1nnnc1SCC(=O)Nc1nc(-c2ccc(Cl)cc2)cs1  1
*n1nnnc1SCc1nc(N)nc(N2CCOCC2)n1 1
*n1s/c(=N/C)nc1-c1ccccc1        1
```

The values of N are much smaller, so the corresponding weight is significantly
smaller.

By default the partition command tries to split the constants evenly (by
weight) across `-n` / `--num-files` files, defaulting to 10, which combined
with the quadratic weighting is why the first few files have only a single,
very common, constant, and why all of the "1" counts are used to fill space in
the remaining files

You can alternatively use `--max-weight` to set an upper bound for the weights
in each file. In this example I'll use the merged fragdb file from the
previous step:

```shell
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.fragdb --max-weight 50000000
Analyzed 'ChEMBL_CYP3A4_hERG.fragdb': #constants: 467865 #fragmentations: 1190215
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG-partition.0000.fragdb' (#1/11, weight: 334589647)
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG-partition.0001.fragdb' (#2/11, weight: 270409141)
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG-partition.0002.fragdb' (#3/11, weight: 225664391)
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG-partition.0003.fragdb' (#4/11, weight: 117895691)
Exporting 10 constants to 'ChEMBL_CYP3A4_hERG-partition.0004.fragdb' (#5/11, weight: 49918518)
Exporting 11 constants to 'ChEMBL_CYP3A4_hERG-partition.0005.fragdb' (#6/11, weight: 49916276)
Exporting 13 constants to 'ChEMBL_CYP3A4_hERG-partition.0006.fragdb' (#7/11, weight: 49899719)
Exporting 7 constants to 'ChEMBL_CYP3A4_hERG-partition.0007.fragdb' (#8/11, weight: 49896681)
Exporting 43 constants to 'ChEMBL_CYP3A4_hERG-partition.0008.fragdb' (#9/11, weight: 49893145)
Exporting 9 constants to 'ChEMBL_CYP3A4_hERG-partition.0009.fragdb' (#10/11, weight: 49879752)
Exporting 467768 constants to 'ChEMBL_CYP3A4_hERG-partition.0010.fragdb' (#11/11, weight: 17615427)
```

If you specify a single fragdb filename then the default output template is
"{prefix}-partition.{i:04}.fragdb" where "{prefix}" is the part of the fragdb
filename before its extension. The idea is to help organize those files
together.

Odds are, you don't want to index the most common fragments. The next two
sections help limits which constants are used.

#### Selecting constants

As you saw, the `mmpdb fragdb_constants` command can be used to list the
constants. It can also be used to list a subset of the constants.

The count for each constant quickly decreases to something a bit more
manageable.

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --limit 20
constant        N
*C      25869
*C.*C   23256
*C.*C.*C        21245
*C.*C.*O        15356
*C.*O   8125
*C.*C.*OC       7076
*C.*OC  6878
*F      6201
*C.*F   6198
*C.*c1ccccc1    5124
*C.*O.*O        5117
*c1ccccc1       5073
*OC     4944
*Cl     4436
*C.*Cl  4388
*O      4300
*F.*F   4281
*C.*F.*F        3935
*C.*C.*F        3656
*F.*F.*F        3496
```

I'll select those constants which occur only 2,000 matches or fewer, and limit
the output to the first 5.

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --max-count 2000 --limit 5
constant        N
*C.*CC.*O       1954
*C.*C(F)(F)F    1915
*C.*C.*OC(C)=O  1895
*C(F)(F)F       1892
*Cl.*Cl 1738
```

or count the number of constants which only occur once (the 1-cut constants
might match with a hydrogen substitution while the others will never match).
I'll use `--no-header` so the number of lines of output matches the number of
constants:

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.fragdb --max-count 1 --no-header | wc -l
  370524
```

These frequent constants are for small fragments. I'll limit the selection to
constants where each part of the constant has at least 5 heavy atoms:

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --min-heavies-per-const-frag 5 --limit 4
constant        N
*c1ccccc1       5073
*c1ccccc1.*c1ccccc1     1116
*Cc1ccccc1      1050
*c1ccc(F)cc1    921
```

I'll also require `N` be between 10 and 1000.

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --min-heavies-per-const-frag 5 \
   --min-count 10 --max-count 1000 --no-header | wc -l
1940
```

That's a much more tractable size for this example.

As you saw earlier, the `mmpdb fragdb_partition` command by default partitions
on all constants. Alternatively, use the `--constants` flag to pass in a list
of constants to use. This can be a file name, or `-` to accept constants from
stdin, as in the following three lines:

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --min-heavies-per-const-frag 5 \
     --min-count 10 --max-count 1000 | \
     mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb --constants -
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG.0000.fragdb' (weight: 423661)
Exporting 1 constants to 'ChEMBL_CYP3A4_hERG.0001.fragdb' (weight: 382376)
Exporting 109 constants to 'ChEMBL_CYP3A4_hERG.0002.fragdb' (weight: 382044)
Exporting 261 constants to 'ChEMBL_CYP3A4_hERG.0003.fragdb' (weight: 382013)
Exporting 261 constants to 'ChEMBL_CYP3A4_hERG.0004.fragdb' (weight: 382013)
Exporting 260 constants to 'ChEMBL_CYP3A4_hERG.0005.fragdb' (weight: 382010)
Exporting 261 constants to 'ChEMBL_CYP3A4_hERG.0006.fragdb' (weight: 382010)
Exporting 262 constants to 'ChEMBL_CYP3A4_hERG.0007.fragdb' (weight: 382010)
Exporting 262 constants to 'ChEMBL_CYP3A4_hERG.0008.fragdb' (weight: 382009)
Exporting 262 constants to 'ChEMBL_CYP3A4_hERG.0009.fragdb' (weight: 382003)
```

Note: the `--constants` parser expects the first line to be a header, which is
why I don't use `--no-header` in the `fragdb_constants` command.
Alternatively, also use `--no-header` in the `fragdb_partition` command if the
input does not have a header.

#### Partitioning in parallel

Partioning large data sets may take significant time because the export
process is single-threaded.

The `fragdb_partition` command can be configured to export only subset of the
partitions using a simple round-robin scheme. If you specify `--task-id n` and
`--num-tasks N` then the given fragdb_partition will only export partitions
`i` such that `i % N == n`.

The expected approach is to create a single constants files which will be
shared by multiple partition commands.

```shell
% mmpdb fragdb_constants ChEMBL_CYP3A4_hERG.*.fragdb --min-heavies-per-const-frag 5 \
     --min-count 10 --max-count 1000 -o constants.dat
```

The following splits the job across two partition commands, with task ids 0
and 1, respectively:

```shell
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb --constants constants.dat --task-id 0 --num-tasks 2
Exporting 1 constants to 'partition.0000.fragdb' (#1/10, weight: 423661)
Exporting 109 constants to 'partition.0002.fragdb' (#3/10, weight: 382044)
Exporting 261 constants to 'partition.0004.fragdb' (#5/10, weight: 382013)
Exporting 261 constants to 'partition.0006.fragdb' (#7/10, weight: 382010)
Exporting 262 constants to 'partition.0008.fragdb' (#9/10, weight: 382009)
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb --constants constants.dat --task-id 1 --num-tasks 2
Exporting 1 constants to 'partition.0001.fragdb' (#2/10, weight: 382376)
Exporting 261 constants to 'partition.0003.fragdb' (#4/10, weight: 382013)
Exporting 260 constants to 'partition.0005.fragdb' (#6/10, weight: 382010)
Exporting 262 constants to 'partition.0007.fragdb' (#8/10, weight: 382010)
Exporting 262 constants to 'partition.0009.fragdb' (#10/10, weight: 382003)
```

Use the `--dry-run` option to get an idea of how many files will be created:
```shell
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.*.fragdb --constants constants.dat --dry-run
i       #constants      weight  filename
0       10      423661  'partition.0000.fragdb'
1       10      382376  'partition.0001.fragdb'
2       10      382044  'partition.0002.fragdb'
3       10      382013  'partition.0003.fragdb'
4       10      382013  'partition.0004.fragdb'
5       10      382010  'partition.0005.fragdb'
6       10      382010  'partition.0006.fragdb'
7       10      382010  'partition.0007.fragdb'
8       10      382009  'partition.0008.fragdb'
9       10      382003  'partition.0009.fragdb'
```
 
#### Indexing in parallel

The partitioned fragdb files can be indexed in parallel:

```shell
% ls partition.*.fragdb | qsub mmpdb index
WARNING: No --output filename specified. Saving to 'partition.0000.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0001.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0002.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0003.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0004.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0005.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0006.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0007.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0008.mmpdb'.
WARNING: No --output filename specified. Saving to 'partition.0009.mmpdb'.
```

(If you don't like these warning messages, use the `--quiet` flag.)

#### Merging partitioned mmpdb files

The last step is to merge the partitioned mmpdb files with the `merge` option,
which only works if no two mmpdb files share the same constant:

```shell
% mmpdb merge partition.*.mmpdb -o ChEMBL_CYP3A4_hERG_distributed.mmpdb
[Stage 1/7] Merging compound records ...
[Stage 1/7] Merged 4428 compound records in 0.046 seconds.
[Stage 2/7] Merging rule_smiles tables ...
[Stage 2/7] Merged 3159 rule_smiles records in 0.030 seconds.
[Stage 3/7] Merging rule tables ...
[Stage 3/7] Merged 21282 rule records in 0.072 seconds.
[Stage 4/7] Merging environment_fingerprint records ...
[Stage 4/7] Merged 1753 environment_fingerprint records in 0.035 seconds.
[Stage 5/7] Merging rule environment records ...
[Stage 5/7] Merged 143661 rule environment records in 0.47 seconds.
[Stage 6/7] Merging constant_smiles and pair records ...
[Stage 6/7] Merged 893 constant SMILES and 203856 pair records in 0.26 seconds
[Stage 7/7] Indexed and analyzed the merged records in 0.33 seconds.
Merged 10 files in 1.3 seconds.
```

Let's take a look:

```shell
% mmpdb list ChEMBL_CYP3A4_hERG_distributed.mmpdb
                Name                 #cmpds #rules #pairs #envs  #stats  |-------- Title --------| Properties
ChEMBL_CYP3A4_hERG_distributed.mmpdb   4428  21282 203856 143661      0  Merged MMPs from 10 files <none>
```

Finally, I'll cross-check this with a normal `mmpdb index`. I need to create
the same subset

```shell
% mmpdb fragdb_partition ChEMBL_CYP3A4_hERG.fragdb --constants constants.dat \
      -n 1 --template ChEMBL_CYP3A4_hERG_subset.fragdb
Exporting 1940 constants to 'ChEMBL_CYP3A4_hERG_subset.fragdb' (#1/1, weight: 3862149)
```

Then index the subset:

```shell
% mmpdb index ChEMBL_CYP3A4_hERG_subset.fragdb
WARNING: No --output filename specified. Saving to 'ChEMBL_CYP3A4_hERG_subset.mmpdb'.
```

And finally, compare the two:

```shell
% mmpdb list ChEMBL_CYP3A4_hERG_subset.mmpdb ChEMBL_CYP3A4_hERG_distributed.mmpdb
                Name                 #cmpds #rules #pairs #envs  #stats  |----------------- Title ------------------| Properties
     ChEMBL_CYP3A4_hERG_subset.mmpdb   4428  21282 203856 143661      0  MMPs from 'ChEMBL_CYP3A4_hERG_subset.fragdb' <none>
ChEMBL_CYP3A4_hERG_distributed.mmpdb   4428  21282 203856 143661      0  Merged MMPs from 10 files                    <none>
```

They are the same, except for the title.


------------------


## History and Acknowledgements


The project started as a fork of the matched molecular pair program
'mmpa' written by Jameed Hussain, then at GlaxoSmithKline Research &
Development Ltd.. Many thanks to them for contributing the code to the
RDKit project under a free software license.

Since then it has gone through two rewrites before the 1.0
release. Major changes to the first version included:

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

The project then forked into three branches:

1. The public GitHub branch, with a few improvements by Christian
  Kramer

2. Andrew Dalke's crowd-funded branch which:
  - replaced the Morgan fingerprint-based hashed environment
  fingerprint with its canonical SMARTS equivalent, and a
  "pseudo-SMILES" which might be used in depictions
  - added Postgres support
  - added export methods to tab-separated and database
  dump formats

3. Mahendra Awale's improvements for:
  - large-database mmpdb generation by partitioning
   on fragment constants
  - playbook generation

Roche funded Andrew Dalke to merge these three branches, resulting in
mmpdb 3.0.

------------------


## Copyright


The mmpdb package is copyright 2015-2023 by F. Hoffmann-La Roche Ltd
and Andrew Dalke Scientific AB, and distributed under the 3-clause BSD
license. See [LICENSE](LICENSE) for details.


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
