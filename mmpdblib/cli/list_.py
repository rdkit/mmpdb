"Implement the 'list' command"

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2019-2021, Andrew Dalke Scientific, AB
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

import click
import sys

from .click_utils import (
    command,
    add_multiple_databases_parameters,
)


list_epilog = """

In the simplest case, look in the current directory for files matching
'*.mmpdb', open each one, and print a terse summary of the information
in the database.

\b
  % mmpdb list
               Name             #cmpds  #rules  #pairs   #envs    #stats   |---------------- Title -----------------| Properties
  CHEMBL_thrombin_Ki_IC50.mmpdb   2985   29513   258294   199657        0  MMPs from 'CHEMBL_thrombin_Ki_IC50.fragdb' <none>
                     csdP.mmpdb   8473 2018581 12372084 12145254 12145254  CSD MP                                     MP

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

\b
  % mmpdb list --all csd.mmpdb
     Name   #cmpds  #rules  #pairs   #envs    #stats   |Title| Properties
  csd.mmpdb   8473 2018581 12372084 12145254 12145254  CSD MP  MP
        Created: 2022-10-08 14:50:16.595104
          #compounds/property:  8473/MP
          #smiles for rules: 56778  for constants: 10759
          Fragment options:
            cut_smarts: [#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]
            max_heavies: 100
            max_rotatable_bonds: 10
            max_up_enumerations: 1000
            method: chiral
            min_heavies_per_const_frag: 0
            min_heavies_total_const_frag: 0
            num_cuts: 3
            rotatable_smarts: [!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]
            salt_remover: <default>
          Index options:
            max_radius: 5
            max_variable_heavies: 10
            min_radius: 0
            smallest_transformation_only: False
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

"""


@command(name="list", epilog=list_epilog)
@click.option(
    "--all",
    "-a",
    "all_option",
    is_flag=True,
    default=False,
    help="List all information about the dataset",
)
@click.option(
    "--recount",
    is_flag=True,
    default=False,
    help="Count the table sizes directly, instead of using cached data",
)
@add_multiple_databases_parameters()
@click.pass_obj
def list_(
    reporter,
    databases_options,
    all_option,
    recount,
):
    """Summarize the contents of zero or more databases

    DATABASES: the mmpdb database files to list (default looks files named '*.mmpdb')
    """

    import json
    from .. import dbutils

    name_list = []
    num_compounds_list = []
    num_rules_list = []
    num_pairs_list = []
    num_envs_list = []
    num_stats_list = []
    titles = []
    property_names = []

    name_width = 5
    num_compounds_width = 6
    num_rules_width = 6
    num_pairs_width = 6
    num_envs_width = 6
    num_stats_width = 6
    title_width = 7

    all_fragment_options = []
    all_index_options = []
    for dbinfo, dataset in dbutils.iter_dbinfo_and_dataset(
        databases_options.databases,
        reporter,
        apsw_warning = False,
    ):
        name = dbinfo.name
        name_width = max(name_width, len(name))
        name_list.append(name)

        table_sizes = dataset.get_table_sizes(recount)
        if not table_sizes.all_defined():
            reporter.warning(
                f"Pre-computed table counts not available in {dbinfo.get_human_name()}. "
                "Forcing --recount."
                )
            table_sizes = dataset.get_table_sizes(recount=True)

        num_compounds = table_sizes.num_compounds
        num_compounds_width = max(num_compounds_width, len(str(num_compounds)))
        num_compounds_list.append(num_compounds)

        num_rules = table_sizes.num_rules
        num_rules_width = max(num_rules_width, len(str(num_rules)))
        num_rules_list.append(num_rules)

        num_pairs = table_sizes.num_pairs
        num_pairs_width = max(num_pairs_width, len(str(num_pairs)))
        num_pairs_list.append(num_pairs)

        num_envs = table_sizes.num_rule_environments
        num_envs_width = max(num_envs_width, len(str(num_envs)))
        num_envs_list.append(num_envs)

        num_stats = table_sizes.num_rule_environment_stats
        num_stats_width = max(num_stats_width, len(str(num_stats)))
        num_stats_list.append(num_stats)

        title = dataset.title
        title_width = max(title_width, len(title))
        titles.append(title)

        prop_names = dataset.get_property_names()
        if prop_names:
            s = " ".join(prop_names)
        else:
            s = "<none>"
        property_names.append(s)

        all_fragment_options.append(dataset.fragment_options_str)
        all_index_options.append(dataset.index_options_str)

    fmt = "%-{}s %-{}s %-{}s %-{}s %-{}s %-{}s  %-{}s Properties".format(
        name_width,
        num_compounds_width,
        num_rules_width,
        num_pairs_width,
        num_envs_width,
        num_stats_width,
        title_width,
    )
    fancy_title = " Title ".center(title_width, "-")
    fancy_title = "|" + fancy_title[1:-1] + "|"
    print(
        fmt
        % (
            "Name".center(name_width),
            "#cmpds".center(num_compounds_width),
            "#rules".center(num_rules_width),
            "#pairs".center(num_pairs_width),
            "#envs".center(num_envs_width),
            "#stats".center(num_stats_width),
            fancy_title,
        )
    )

    fmt = "%{}s %{}d %{}d %{}d %{}d %{}d  %-{}s %s".format(
        name_width,
        num_compounds_width,
        num_rules_width,
        num_pairs_width,
        num_envs_width,
        num_stats_width,
        title_width,
    )

    prefix = " " * num_compounds_width
    for (
        name,
        num_compounds,
        num_rules,
        num_pairs,
        num_envs,
        num_stats,
        title,
        names_and_counts,
        fragment_options,
        index_options,
    ) in zip(
        name_list,
        num_compounds_list,
        num_rules_list,
        num_pairs_list,
        num_envs_list,
        num_stats_list,
        titles,
        property_names,
        all_fragment_options,
        all_index_options,
    ):
        print(
            fmt
            % (
                name,
                num_compounds,
                num_rules,
                num_pairs,
                num_envs,
                num_stats,
                title,
                names_and_counts,
            )
        )
        if all_option:
            creation_date = dataset.creation_date
            creation_date_str = creation_date.isoformat(" ")
            print(prefix + "Created:", creation_date_str)

            s = " "  # Always have a trailing space
            for property_name, count in dataset.get_property_names_and_counts():
                s += "%s/%s " % (count, property_name)
            if s == " ":
                s = "(no properties)"
            else:
                s = s[:-1]  # strip the trailing space
            print(prefix + "  #compounds/property:", s)

            print(
                prefix
                + "  #smiles for rules: %d  for constants: %d"
                % (dataset.get_num_rule_smiles(), dataset.get_num_constant_smiles())
            )

            options = json.loads(fragment_options)
            print(prefix + "  Fragment options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))

            options = json.loads(index_options)
            print(prefix + "  Index options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))
