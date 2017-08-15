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

import argparse

from . import __version__ as _version
from . import config
from . import smarts_aliases
from . import do_help


#### Help construct the argparser

from .config import nonnegative_int, cutoff_list

def add_single_dataset_arguments(parser):
    parser.add_argument("databases", nargs=1, metavar="DATABASE",
                        help="location of the database (eg, 'csd.mmpdb')")
    
def add_multiple_dataset_arguments(parser):
    parser.add_argument("databases", nargs="*", metavar="DATABASE",
                        help="zero or more database locations (eg, 'csd.mmpdb')")

def add_single_property_arguments(parser):
    parser.add_argument("--property", "-p", metavar="NAME", required=True,
                        help="property to use")

def add_multiple_properties_default_all(parser):
    parser.add_argument("--property", "-p", metavar="NAME", action="append",
                        help="property to use")
    parser.add_argument("--no-properties", action="store_true",
                        help="don't use any properties")
    parser.set_defaults(all_properties=False)

def add_in_memory(parser):
    p.add_argument("--in-memory", action="store_true",
                    help="load the SQLite database into memory before use (requires APSW)")

####

parser = argparse.ArgumentParser(
    description="Matched-molecular pair database loader",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
The 'mmpdb' program implements a set of subcommands to work with
a matched molecular pair database. The commands are roughly
categorized as "analysis" and "admin" commands.

The analysis commands fragment a SMILES database, index the fragments
into matched molecular pairs into a local SQLite database, import
molecular properties into a database, and searches the database for
possible transformations or predictions.

The admin commands are used to administer a database. They include
include ways to list the available data sets, dump the data as
a SMILES file or CSV file, update new properties, and
re-aggregate the rule dataset should any property values change.

For a short description of how to generate and use a dataset:
  % mmpdb help-analysis

For a short description of how to adminster a database:
  % mmpdb help-admin

The "help-*-format" commands (like "help-property-format") give
more details about the given format.

In addition, pass the "--help" option to a given command to see
the full list of options for the command.

""")

parser.add_argument("--quiet", "-q", action="store_true",
                    help="do not show progress or status information")
parser.add_argument('--version', action="version", version="%(prog)s " + _version)

subparsers = parser.add_subparsers()


#### mmpdb fragment

p = fragment_parser = subparsers.add_parser(
    "fragment",
    help="fragment structures in a SMILES file based on its rotatable bonds",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Fragment molecules in a SMILES file by breaking on 'cut bonds',
as matched by --cut-smarts. Cut up to --num-cuts bonds. Don't
fragment molecules with more than --max-rotatable-bonds bonds or
--max-heavies heavy atoms.

The input structures comes from a SMILES file. By default the fields
are whitespace delimited, where the first column contains the SMILES
and the second contains the identifier. Use --delimiter to change
the delimiter type, and --has-header to skip the first line. See
"mmpdb help-smiles-format" for more details.

The input SMILES strings are pre-processed to remove salts before
fragmenting. The default uses the default RDKit SmilesRemover. Use
--salt-remover to specify alternative rules, or use the special
name '<none>' to not remove any salts.

By default the fragmentation method uses 4 threads, which gives a
nearly 4-fold speedup. Use --num-jobs change the number of threads.

It can take a while to generate fragments. Suppose you want to update
the compound set every week, where only a few records are added or
changed. Most of the fragments will be the same between those two data
sets. What you can do is specify the old fragments file as a --cache so
the fragmentation method can re-use the old fragmentation, assuming
the structure hasn't changed for a given record.

Examples:

1) Fragment the SMILES file to produce a fragments file

  % mmpdb fragment CHEMBL_thrombin_Ki_IC50.smi -o CHEMBL_thrombin_Ki_IC50.fragments

2) Read from a gzip-compressed tab-delimited SMILES file. Use 8
threads to fragment the structures. Save the results to
dataset.fragments.gz .

  % mmpa fragment --delimiter tab dataset.smi.gz --num-jobs 8 \\
      -o dataset.fragments.gz

3) Fragment the SMILES in 'dataset.smi.gz'. Reuse fragment information
from the cache file 'old_dataset.fragments.gz' if possible, instead of
computing the fragments from scratch each time. Save the results to
'new_dataset.fragments.gz'.

  % mmpa fragment --cache old_dataset.fragments.gz dataset.smi.gz \\
      -o new_dataset.fragments.gz


""" + smarts_aliases.get_epilog("--cut-smarts", smarts_aliases.cut_smarts_aliases)
)



def fragment_command(parser, args):
    from . import do_fragment
    do_fragment.fragment_command(parser, args)

config.add_fragment_arguments(p)

p.add_argument("--cache", metavar="SOURCE",
               help="get fragment parameters and previous fragment information from SOURCE")
p.add_argument("--num-jobs", "-j", metavar="N", type=config.positive_int, default=4,
               help="number of jobs to process in parallel (default: 4)")

p.add_argument("-i", "--in", metavar="FORMAT", dest="format",
               choices=("smi", "smi.gz"),
               help="input structuture format (one of 'smi', 'smi.gz')")
p.add_argument("--delimiter", default="whitespace",
               # 'native' is hidden support for for chemfp compatability
               choices=("whitespace", "to-eol", "comma", "tab", "space", "native"),
               help="SMILES file delimiter style (one of 'whitespace' (default), 'to-eol', 'comma', 'tab', or 'space')")
p.add_argument("--has-header", default=False, action="store_true",
               help="skip the first line, which is the header line")
p.add_argument("--output", "-o", metavar="FILENAME",
               help="save the fragment data to FILENAME (default=stdout)")
p.add_argument("--out", metavar="FORMAT", choices=("fragments", "fragments.gz", "fraginfo", "fraginfo.gz"),
               help="output format. One of 'fragments' or 'fragments.gz'. "
               "If not present, guess from the filename, and default to 'fragments'")
p.add_argument("structure_filename", nargs="?", default=None,
               help="SMILES filename (default: read from stdin)")

    
p.set_defaults(command=fragment_command,
               subparser=p)

#### mmpdb smifrag

p = smifrag_parser = subparsers.add_parser(
    "smifrag",
    help="fragment a SMILES and print the each variable and constant fragment",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Fragment a SMILES and print details about each variable and constant
fragment and how they are connected.

""" + smarts_aliases.get_epilog("--cut-smarts", smarts_aliases.cut_smarts_aliases)
)

def smifrag_command(parser, args):
    from . import do_fragment
    do_fragment.smifrag_command(parser, args)


config.add_fragment_arguments(p)
p.set_defaults(command=smifrag_command,
               subparser=p)
p.add_argument("smiles", metavar="SMILES",
               help="SMILES string to fragment")

#### mmpdb index

p = index_parser = subparsers.add_parser(
    "index",
    help="index fragments and find matched molecular pairs",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Read a fragments file, match the molecular pairs, and save the results
to a mmpdb database file. Use "mmpdb help-analysis" for an explanation
about the terminology and process. A matched molecular pair connects
structure 1, made of variable part 1 (V1) and constant part C, with
structure 2, made of variable part 2 (V2) and constant part C, via the
transformation V1>>V2. The following uses |X| to mean the number of
heavy (non-isotopic hydrogen) atoms in X.

There are several ways to restrict which fragmentations or pairings
will to allow. The --min-variable-heavies and --max-variable-heavies
options set limits to |V1| and |V2|. The --min-variable-ratio and
--max-variable-ratio set limits on |V1|/|V1+C| and |V2|/|V2+C|. The
--max-heavies-transf option sets a maximum bound on abs(|V1|-|V2|).

The --max-frac-trans places a limit on the number of atoms which take
part in the transformation, defined as |V+C(r)|/|V+C| where C(r) is
the number of atoms in the circular environment of the constant part
for a given radius r. C(0) is 0. The goal is to be able to exchange
transforms and environment fingerprints but minimize the relative
amount of information revealed about the constant part.

The filter --max-variable-heavies is always used, with a default value
of 10. Specify "none" if you want no limit.

The transformation can be described as V1>>V2 or V2>>V1. By default
only one transformation is output; the one with the smallest value,
alphabetically. The --symmetric flag outputs both directions.

The optional --properties file contains physical property information
in a table format. Use "mmpdb help-property-format" for details. If a
property file is specified during the indexing stage then only those
compounds with at least one property are indexed. See "mmpdb
loadprops" for a way to specify properties after indexing.

By default the output will be saved in 'mmpdb' format using a filename
based on the input fragment filename and with the extension
'.mmpdb'. Use --output to specify an output filename. Use --out to
change the output format. By default the alternate formats will write
to stdout.

The 'mmpdb' format is based on a SQLite database. The 'mmpa' format
stores the same data in a text file with tab-separated fields. The
'csv' format is a tab-separated table with the columns:
  SMILES1  SMILES2  id1  id2  V1>>V2  C

The csv format does not include property information. The mmpa and csv
format also support gzip compression.

The --title specifies a string which goes into the database. The idea
is to store a label or short description to be displayed to the
user. If not given, it uses a short description based on the input
filename.

The --memory option writes a summary of memory use to stderr if the
'psutil' package is available. (If not, it prints a warning message
that it needs the package.) This was used during development to figure
out ways to reduce overall memory use.

The experimental --in-memory option loads the SQLite database into
memory before using it. This option requires the third-party APSW
SQLite adapater and will not work with Python's built-in SQLite
module. The transformation analysis does many scattered database
lookups. In some cases, like with a network file system, it can be
faster to load the entire database into memory where random access is
fast, rather than do a lot of disk seeks.


Examples:

1) Index 'csd.fragments' and save the results to 'csd.mmpdb'. The
title will be "MMPs from 'csd.fragments'".

  % mmpdb index csd.fragments

2) Index 'csd.fragments', use properties from 'MP.csv' and limit the
pair matching to compounds listed in the CSV file, set the title to
'CSD MP', and save the results to 'csd_MP.mmpdb'.

  % mmpdb index csd.fragments --properties MP.csv --title "CSD MP" -o csd_MP.mmpdb

3) Limit the indexing to variable terms which have at least 12 heavy
atoms and where the size of the variable is no more than 40% of the
entire structure, and save the transformation in both A>>B and B>>A
(symmetric) forms:

  % mmpdb index CHEMBL_Thrombin_logD.fragments --symmetric \\
      --max-variable-ratio 0.4 --max-variable-heavies 12 \\
      --title "CHEMBL ratio 40%" --output CHEMBL_ratio_40.mmpdb
""")

def index_command(parser, args):
    from . import do_index
    do_index.index_command(parser, args)

config.add_index_options(p)
p.add_argument("--symmetric", "-s", action="store_true",
               help="Output symmetrically equivalent MMPs, i.e output both cmpd1,cmpd2, "
                      "SMIRKS:A>>B and cmpd2,cmpd1, SMIRKS:B>>A")
p.add_argument("--properties", metavar="FILENAME",
               help="File containing the identifiers to use and optional physical properties")
p.add_argument("--output", "-o", metavar="FILENAME",
               help=("save the fragment data to FILENAME. "
                     "Default for mmpdb is based on the fragment filename, "
                     "otherwise stdout."))
p.add_argument("--out", metavar="FORMAT", choices=("csv", "csv.gz", "mmpa", "mmpa.gz", "mmpdb"),
               help="Output format. One of 'mmpdb' (default), 'csv', 'csv.gz', 'mmpa' or 'mmpa.gz'. "
               "If not present, guess from the filename, and default to 'mmpdb'")
p.add_argument("--title",
               help="A short description of the dataset. If not given, base the title on the filename")
p.add_argument("--memory", action="store_true",
               help="Report a summary of the memory use")
p.add_argument("fragment_filename", nargs="?", default=None,
               help="SMILES filename (default: read from stdin)")

p.set_defaults(command=index_command,
               subparser=p)

######## Work with a database


#### mmpdb list
p = list_parser = subparsers.add_parser(
    "list",
    help="summarize the contents of zero or more databases",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

In the simplest case, look in the current directory for files matching
'*.mmpdb', open each one, and print a terse summary of the information
in the database.

  % mmpdb list
               Name             #cmpds  #rules  #pairs   #envs    #stats   |-------- Title ---------| Properties
  CHEMBL_thrombin_Ki_IC50.mmpdb   2985   29513   258294   199631        0  thrombin Ki IC50           
                      csd.mmpdb  16843 2525166 15328608 15199376 15199376  MMPs from 'csd.fragments'  MP

The output is a set of columns. The first line is the header. The first
column contains the database name. The next columns contain the number
of compounds, number of rules, number of pairs (a rule may have many
matched molecular pairs), number of rule environments (a rule may have
many environments), and number of property statistics for the rule
environments. After that is the user-defined title field, followed by
a list of the property or activity names stored.

The first entry, for thrombin, has no properties, which is why it also
has no property statistics. The second entry has a 'MP' property,
which in this case means 'melting point'.

The specific database location(s) can be given on the
command-line. The '--all' option shows more detailed information about
the dataset. The following gives more detailed information about the
database 'csd.mmpdb':


  % mmpdb list --all csd.mmpdb 
     Name   #cmpds  #rules  #pairs   #envs    #stats   |-------- Title --------| Properties
  csd.mmpdb  16843 2525166 15328608 15199376 15199376  MMPs from 'csd.fragments' MP
        Created: 2017-05-26 15:07:07.775249
          #compounds/property:  16843/MP
          #smiles for rules: 81538  for constants: 21520
          Fragment options:
            cut_smarts: [#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]
            max_heavies: 100
            max_rotatable_bonds: 10
            method: chiral
            num_cuts: 3
            rotatable_smarts: [!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]
            salt_remover: <default>
          Index options:
            max_variable_heavies: 10
            symmetric: False

'Created' shows the creation time. '#compounds/property' shows how
many compounds have a given property, for each of the available
properties. The '#smiles' line says how many distinct SMILES strings
are used for the rules and the constants tables. 'Fragment options'
and 'Index options' are, I think, self-explanatory.

The count fields (like the number of compounds and rules) are
pre-computed and stored in the database. If the database is updated
incorrectly, it is possible for the cached information to be
invalid. Use '--recount' to have SQLite compute the values directly
from the database contents.

""")

def list_command(parser, args):
    from . import do_database
    do_database.list_command(parser, args)


p.add_argument("--all", "-a", action="store_true",
               help="list all information about the dataset")
p.add_argument("--quiet", "-q", action="store_true",
               help="do not show progress or status information")
p.add_argument("--recount", action="store_true",
               help="count the table sizes directly, instead of using cached data")
add_multiple_dataset_arguments(p)
p.set_defaults(command=list_command,
               subparser=p)




#### mmpdb loadprops

p = loadprops_parser = subparsers.add_parser(
    "loadprops",
    help="load properties for existing structures",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Load structure property values from a CSV into a data set.

The property file contains a data table. Here is an example:

  ID MP CHR1 CHR2
  GEJYOJ 3 71 31.3
  ACIDUL 5 65 67.2
  KIXRIS 5 * *
  SOFWIV01 5 83 79.3

The fields should be tab separated. For full details use
"mmpdb help-property-format".

If the given property value already exists in the database then the
existing database value will be updated. Otherwise loadprops will
create a new property record. If an identifier isn't in the database
then its values will be ignored.

After importing the data, the corresponding aggregate values for the rules
will be recalculated.

Example:

  % mmpdb loadprops --properties MP.csv mmpdb.db
  
"""
)

def loadprops_command(parser, args):
    from . import do_database
    do_database.loadprops_command(parser, args)

add_single_dataset_arguments(p)
p.add_argument("--properties", "-p", nargs="?", metavar="FILENAME",
               help="File containing the identifiers to use and optional physical properties")
p.set_defaults(command=loadprops_command,
               subparser=p)
    

#### mmpdb reaggregate

# Commented out because it doesn't seem useful
'''
p = reaggregate_parser = subparsers.add_parser(
    "reaggregate",
    help="recompute the rule statistics",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Compute the rule statistics using the property values in the
database. Normally this will not change anything because 'mmpdb index'
(if --properties are given) and 'mmpdb loadprops' will compute the
statistics automatically, and those are the only two commands which
can change rule properties.

This might be useful if you update property values in the database
yourself, and want mmpdb to handle the update of the aggregate
statistics.

""")

def reaggregate_command(parser, args):
    from . import do_database
    do_database.reaggregate_command(parser, args)

add_single_dataset_arguments(p)
p.add_argument("--property", "-p", metavar="NAME", action="append",
               help="property to use (default reaggregates all properties)")
p.set_defaults(command=reaggregate_command,
               subparser=p, no_properties=False, all_properties=False)
'''

#### mmpdb transform


p = transform_parser = subparsers.add_parser(
    "transform",
    help="transform a structure",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Apply transforms from an mmpdb database to an input structure.
Include possible property change statistics for the resulting
products.

Specify the input structure using --smiles. This will be fragmented
using the fragmentation parameters appropriate for the database. By
default all fragmentations will be considered. Use --min-variable-size
and --min-constant-size to set minimum heavy counts for the constant
and variable parts of the fragment.

By default the matching algorithm evaluates all radii around the local
environment of the constant's connection points. The scoring function
(see below) decides which radius is best. Use the --min-radius option
to require the environment match up to at least N bonds away.

Use --min-pairs to require that a transformation have at least N pairs
in the database. The default is 0, which allows all transformations.

For more complex filters, use the --where option. It takes a Python
expression which is allowed to use any of the following variable
names:

  rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles,
  to_num_heavies, smirks, rule_environment_id, radius, fingerprint_id,
  fingerprint, rule_environment_statistics_id, count, avg, std,
  kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value,
  is_bidirectional

as well as the Python variables None, True, and False. The values come
from the transformation rule, the rule environment, and the rule
environment statistics.

The --where expression is evaluated independently for each
property. There is currently no way to specify a different expression
for different properties, or for the expression to know which property
is being evaluated.

Sometimes multiple transformations of an input structure lead to the
same final product. When this happens, the possible transformation are
put into bins based on their number of pairs. By default the first bin
contains transforms with at least 10 pairs, the second contains
transforms with at least 5 pairs, and the third contains transforms
with any pairs. The bin thresholds can be changed with the
--rule-selection-cutoffs option, which contains a list of
comma-separated integers. The default is equivalent to
'--rule-selection-cutoffs 10,5,0'.

A scoring function is used to decide which transformation to use from
a bin. The transformation with the largest score is selected. Use the
--score option to define an alternate scoring function as a Python
expression. The default is equivalant to:

  ((ninf if std is None else -std), radius, from_num_heavies, from_smiles)

This expression may contain any of the variables in the --where
option, as well as 'inf' and 'ninf', which can be used as equivalent
for positive and negative infinity in numerical comparisons. The
default expression selects the smallest standard deviation, breaks
ties using the largest radius, breaks ties with the most atoms
changed, then breaks ties arbitrarily by the substructure SMILES of
the left side of the transformation.

Note: if the --where or --score expressions start with a '-' then the
command-line parser may confuse it with a command-line option. In that
case, use a space as the first character in the expression, or enclose
the expression with parentheses.

Specify a SMARTS pattern with --substructure to limit the output to
products containing the specified substructure.

By default the output will contain predicted property changes for all
of the properties in the database. Use --property to select specific
options. Use it once for each property you want to evaluate. For
example, to get the results for both MW and MP use:
  --property MW --property MP

If you do not want property information in the output, specify
--no-properties.

The transformation output by default is sent to stdout. Use '--output'
to specify an output filename.

The output is a tab-delimited CSV file where the first line is a
header. The first column contains a sequential identifier, and the
second column contains the SMILES string for the product. Next come
the property columns. Each property gets 17 columns of output. The
column headers are prefixed with the property name followed by a "_"
and then the name of the data stored in the column. (There is
currently no way to limit the output to specific columns, other than
to change the code in do_transform.py.)

The "*_from_smiles" and "*_to_smiles" columns describe the
transformation.  Sometimes multiple transforms lead to the same
product and cause different properties to have different
transformations. The "*_radius" and "*_fingerprint" columns contain
the environment fingerprint radius and fingerprint string. The
remaining columns contain the database row id for the rule environment
record, and the statistics information for the given property.

The '--explain' option writes debug information to stderr. The
'--times' options reports timing information about the major stages to
stderr.

The transform code runs in a single thread by default. It has been
parallelized, though it is not highly scalable. Use '--jobs' to
specify the number of threads (really, "processes") to use.

Examples:

1) Generate all of the products of diphenyl ether using the MMP
transforms in the 'csd.mmpdb' database where there are are least 30
pairs. Also include the predicted effects on the 'MP' property. (Note:
the output is reformatted and trimmed for use as help text.)

  % mmpdb transform csd.mmpdb --smiles 'c1ccccc1Oc1ccccc1' --min-pairs 30 -p MP
  ID                 SMILES MP_from_smiles       MP_to_smiles  MP_radius  \\
   1  Brc1ccc(Oc2ccccc2)cc1  [*:1]c1ccccc1  [*:1]c1ccc(Br)cc1          0
   2  COc1ccc(Oc2ccccc2)cc1  [*:1]c1ccccc1  [*:1]c1ccc(OC)cc1          0
   3             COc1ccccc1  [*:1]c1ccccc1             [*:1]C          0

                               MP_fingerprint  MP_rule_environment_id  \\
  59SlQURkWt98BOD1VlKTGRkiqFDbG6JVkeTJ3ex3bOA                     947
  59SlQURkWt98BOD1VlKTGRkiqFDbG6JVkeTJ3ex3bOA                    4560
  59SlQURkWt98BOD1VlKTGRkiqFDbG6JVkeTJ3ex3bOA                      90

  MP_count   MP_avg  MP_std  MP_kurtosis  MP_skewness  MP_min  MP_q1  \\
        34  14.5290  30.990    -0.267780      0.32663     -66   -7.0
        56   8.7143  38.945     7.013600      1.81870    -172  -10.0
       106 -23.4430  36.987     1.563800      0.65077    -159  -44.0

  MP_median  MP_q3  MP_max  MP_paired_t    MP_p_value
       15.5   37.0      67      -2.7338  9.987200e-03
       10.5   32.5      79      -1.6745  9.971500e-02
      -20.0   -3.0      49       6.5256  2.447100e-09

2) Require a standard deviation of no larger than 4.5 and give
priority to transformation with at least 20 pairs before following the
normal cutoffs. The --score here matches the default scoring function.

  % mmpdb transform csd.mmpdb --smiles 'c1ccccc1Oc1ccccc1' \\
      --where 'std is not None and std < 4.5' \\
      --rule-selection-cutoffs '20,10,5,0' \\
      --score '((ninf if std is None else -std), radius, from_num_heavies, from_smiles)' \\
      --property MP


""")

def transform_command(parser, args):
    from . import do_analysis
    do_analysis.transform_command(parser, args)

add_single_dataset_arguments(p)
p.add_argument("--smiles", "-s", required=True,
               help="the base structure to transform")
p.add_argument("--min-variable-size", type=nonnegative_int, metavar="N", default=0,
               help="require at least N atoms in the variable fragment (default: 0)")
p.add_argument("--min-constant-size", type=nonnegative_int, metavar="N", default=0,
               help="require at least N atoms in the constant fragment (default: 0)")
p.add_argument("--min-radius", "-r", choices=("0", "1", "2", "3", "4", "5"), default="0",
               help="fingerprint radius (default: 0)")
p.add_argument("--min-pairs", type=nonnegative_int, metavar="N", default=0,
               help="require at least N pairs in the transformation to report a product (default: 0)")
p.add_argument("--substructure", "-S", metavar="SMARTS",
               help="require the substructure pattern in the product")
add_multiple_properties_default_all(p)
config.add_rule_selection_arguments(p)
add_in_memory(p)
p.add_argument("--jobs", "-j", type=config.positive_int, default=1,
               help="number of jobs to run in parallel (default: 1)")
p.add_argument("--explain", action="store_true",
               help="explain each of the steps in the transformation process")
p.add_argument("--output", "-o", metavar="FILENAME",
               help="save the output to FILENAME (default=stdout)")
p.add_argument("--times", action="store_true",
               help="report timing information for each step")
p.set_defaults(command=transform_command,
               subparser=p)

    
#### mmpdb predict

p = predict_parser = subparsers.add_parser(
    "predict",
    help="predict the effect of a structural transformation",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

  --smiles SMILES, -s SMILES
                        the base structure to transform
  --reference SMILES    the reference structure
  --property NAME, -p NAME
                        property to use

  --value VALUE, -v VALUE
                        the property value for the reference

Predict the change in a specified '--property' going from the
'--reference' SMILES to the input '--smiles', using the matched
molecular pair 'DATABASE.' If several rules lead to the same input
'--smiles', use the one prioritized by '--rule-selection-cutoffs' and
'--score'.

By default the predicted delta value and standard deviation will be
written to stdout, in the format:
  predicted delta: -4.91667 +/- 18.0477

If the optional reference property '--value' is specified then the
predicted output value will also be included, in the format: predicted
delta: -4.91667 predicted value: 7.08333 +/- 18.0477

The '--where' option takes a Python expression which describes which
transforms to consider. The '--rule-selection-cutoffs' and '--score'
options select which transform to use if multiple transforms produce
the same product. If the third-party APSW module is available then the
'--in-memory' option loads the SQLite database into memory before
doing the prediction. The '--explain' and '--times' options send
respectively debug and timing information to stderr. More details are
available from 'mmpdb transform --help'.

Use '--save-details' to save details about the possible transformations
to two tab-separated CSV files. The filenames are of the form
${PREFIX}_rules.txt and ${PREFIX}_pairs.txt, where ${PREFIX} is
specified by --prefix and by default is 'pred_details'.

The '${PREFIX}_rules.txt' file contains the following columns for each
transform, based on the rule, rule environment, and rule environment
statistics:
  rule_id, rule_environment_id, radius, fingerprint, from_smiles,
  to_smiles, count, avg, std, kurtosis, skewness, min, q1, median, q3,
  max, paired_t, p_value

The '${PREFIX}_pairs.txt' file contains the following columns for
each of the pairs in each of the transforms:
  from_smiles, to_smiles, radius, fingerprint, lhs_public_id,
  rhs_public_id, lhs_smiles, rhs_smiles, lhs_value, rhs_value, delta


Examples:

1) Predict the effect of substituting a sulfur in diphenyl ether if
the known melting point is 12C:

  % mmpdb predict csd.mmpdb --smiles c1ccccc1Sc1ccccc1 \\
            --reference c1ccccc1Oc1ccccc1 --property MP --value 12.0
  predicted delta: +4.91667 predicted value: 16.9167 +/- 18.0477

2) Do the same calculation and save details about the transformations
to O_to_S_rules.txt and O_to_S_pairs.txt:

  % mmpdb predict csd.mmpdb --smiles c1ccccc1Sc1ccccc1 \\
            --reference c1ccccc1Oc1ccccc1 --property MP --value 12.0 \\
            --save-details --prefix O_to_S

"""
    )
# predict a property for a structure given the known property for
# another structure, using the MMPA database to identify the possible
# transformations and predicted property changes. The input to the
# tool will be the reference structure, the prediction structure, and
# the properties to predict. There will also be an option to report
# pairs with the same transformations and their properties.

def predict_command(parser, args):
    from . import do_analysis
    do_analysis.predict_command(parser, args)

add_single_dataset_arguments(p)
p.add_argument("--smiles", "-s", metavar="SMILES", required=True,
               help="the base structure to transform")
p.add_argument("--reference", metavar="SMILES", required=True,
               help="the reference structure")
add_single_property_arguments(p)
config.add_rule_selection_arguments(p)
p.add_argument("--value", "-v", type=float, default=None,
               help="the property value for the reference")
add_in_memory(p)
p.add_argument("--explain", action="store_true",
               help="explain each of the steps in the prediction process")
p.add_argument("--save-details", action="store_true",
               help="save information about the transformation pairs and statistics to two CSV files")
p.add_argument("--prefix", metavar="STRING", default="pred_detail",
               help="prefix to use for each CSV filename (default: 'pred_details')")
p.add_argument("--times", action="store_true",
               help="report timing information for each step")

p.set_defaults(command=predict_command,
               subparser=p)


##### mmpdb smicat

p = smicat_parser = subparsers.add_parser(
    "smicat",
    help="write the database compounds to a SMILES file",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Write information about the compounds in DATABASE, formatted as a
SMILES file.

Each compound has two associated SMILES, the input SMILES used as
input to fragmentation, and the canonical SMILES string from RDKit
after input processing (typically desalting and structure
normalization). By default the output uses the processed SMILES. Use
'--input-smiles' to use the input SMILES string.

By default the output SMILES file is written to stdout. Use '--output'
to save the output to the named file.

Examples:

1) Write the cleaned-up structures as a SMILES file to stdout:

  % mmpdb smicat csd.mmpdb

2) Save the structures to the file "original.smi", and use the input
SMILES instead of the de-salted SMILES:

  % python mmpdb smicat csd.mmpdb -o original.smi --input-smiles

""")

def smicat_command(parser, args):
    from . import do_database
    do_database.smicat_command(parser, args)

add_single_dataset_arguments(p)
p.add_argument("--input-smiles", action="store_true",
               help="Use the input SMILES instead of the cleaned-up SMILES")
p.add_argument("--output", "-o",
               help="output filename (default is stdout)")
p.set_defaults(command=smicat_command,
               subparser=p)

##### mmpdb propcat

p = propcat_parser = subparsers.add_parser(
    "propcat",
    help="write the database properties to a properties file",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Write information about the properties for the compounds in DATABASE,
formatted as a property file. Use `mmpdb help-property-file` for
details about the property file format.

The output from this command is a tab-delimited CSV file where the
first column has the head "ID" and contains the compound identifier.
The other columns contain property information for each compound. The
column title is the property name.

By default there is one column for each property in the databases, and
the one row for each compound with at least one property. Use
'--property' to limit the output to a specific property, or use it
multiple times to specify multiple property names to output. Use
'--all' to list all of the compounds, even if the compound has none of
the specified properties.

The character "*" will be use if a listed compound is missing a given
property.

Examples:

1) Write all of the properties to stdout:

  % mmpdb propcat CHEMBL_thrombin_Ki_IC50.mmpdb

2) Write the "MP" property to "MP.properties":

  % mmpdb propcat csd.mmpdb --property MP -o MP.properties

3) Write the compound identifiers only to stdout:

  % mmpdb propcat csd.mmpdb --no-properties --all

""")

def propcat_command(parser, args):
    from . import do_database
    do_database.propcat_command(parser, args)

add_single_dataset_arguments(p)
add_multiple_properties_default_all(p)
p.add_argument("--all", action="store_true",
               help="include compounds which have no properties")
p.add_argument("--output", "-o",
               help="output filename (default is stdout)")
p.set_defaults(command=propcat_command,
               subparser=p)

#### mmpdb drop_index
p = drop_index_parser = subparsers.add_parser(
    "drop_index",
    help="drop the database indices",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Drop the database indices from DATABASE. This is mostly used during
development. The index takes about 1/2 of the size of the database, so
if you need to save space for data exchange or archival purposes then
you might drop the indices, and re-create them later when needed.
""")

def drop_index_command(parser, args):
    from . import do_database
    do_database.drop_index_command(parser, args)

add_single_dataset_arguments(p)
p.set_defaults(command=drop_index_command,
               subparser=p)


#### mmpdb create_index
p = create_index_parser = subparsers.add_parser(
    "create_index",
    help="create the database indices",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""

Create the database indices for DATABASE. This is mostly used during
development.
""")

def create_index_command(parser, args):
    from . import do_database
    do_database.create_index_command(parser, args)

add_single_dataset_arguments(p)
p.set_defaults(command=create_index_command,
               subparser=p)



##### Help

do_help.add_help_commands(subparsers)

def main(argv=None):
    # if there are extra arguments in the subparser then
    #   args = parser.parse_args(argv)
    # will pass the extra fields back to the main parser,
    # which generates an un-informative top-level error
    # message instead of a more specific subparser error.
    # Instead, I get the known arguments, and send any 
    # remaining subarguments to the correct subparser.
    # (See http://stackoverflow.com/questions/25333847/argparse-subcommand-error-message
    # for more details.)
    parsed_args, remaining_argv = parser.parse_known_args(argv)
    if remaining_argv:
        subparser = getattr(parsed_args, "subparser", parser)
        subparser.error("unrecognized arguments: %s" % (" ".join(remaining_argv)))

    # Python 2 raises an exception. Python 3 doesn't.
    # This means I can provide a slightly better message.
    command = getattr(parsed_args, "command", None)
    if command is None:
        parser.error("too few arguments. Use --help for more information.")
        
        
    parsed_args.command(parsed_args.subparser, parsed_args)

if __name__ == "__main__":
    main()
