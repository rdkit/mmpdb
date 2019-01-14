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

# list [database*]

import sys
import json
import itertools
#import time

from . import schema
from . import command_support
from . import dbutils
from . import index_algorithm
from . import fileio    

# mmpdb list [--all] [--quiet] [--recount] [filename]*

def list_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)
    databases = args.databases

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
    for dbinfo in dbutils.iter_dbinfo(databases, reporter):
        reporter.update("Opening %r ... " % (dbinfo.get_human_name(),))
        database = None
        try:
            database = dbinfo.open_database()
            dataset = database.get_dataset()
        except dbutils.DBError as err:
            reporter.update("")
            reporter.report("Skipping %s: %s"
                            % (dbinfo.get_human_name(), err))
            if database is not None:
                database.close()
            continue
        reporter.update("")

        name = dbinfo.name
        name_width = max(name_width, len(name))
        name_list.append(name)

        table_sizes = dataset.get_table_sizes(args.recount)
        
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
        name_width, num_compounds_width, num_rules_width, num_pairs_width,
        num_envs_width, num_stats_width, title_width)
    fancy_title = " Title ".center(title_width, "-")
    fancy_title = "|" + fancy_title[1:-1] + "|"
    print(fmt % ("Name".center(name_width), "#cmpds".center(num_compounds_width), "#rules".center(num_rules_width),
                 "#pairs".center(num_pairs_width), "#envs".center(num_envs_width), "#stats".center(num_stats_width),
                 fancy_title))
    
    fmt = "%{}s %{}d %{}d %{}d %{}d %{}d  %-{}s %s".format(
        name_width, num_compounds_width, num_rules_width, num_pairs_width, num_envs_width, num_stats_width, title_width)
    
    prefix = " "*num_compounds_width
    for (name, num_compounds, num_rules, num_pairs, num_envs, num_stats, title,
         names_and_counts, fragment_options, index_options) in zip(
            name_list, num_compounds_list, num_rules_list, num_pairs_list, num_envs_list, num_stats_list, titles,
            property_names, all_fragment_options, all_index_options):
        print(fmt % (name, num_compounds, num_rules, num_pairs, num_envs, num_stats, title, names_and_counts))
        if args.all:
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
            
            print(prefix + "  #smiles for rules: %d  for constants: %d"
                  % (dataset.get_num_rule_smiles(), dataset.get_num_constant_smiles()))
            
            options = json.loads(fragment_options)
            print(prefix + "  Fragment options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))
    
            options = json.loads(index_options)
            print(prefix + "  Index options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))
    

# mmpdb create_index <filename>
def create_index_command(parser, args):
    mmpdb = dbutils.open_database_from_args_or_exit(args)
    with mmpdb.atomic():
        schema.create_index(mmpdb)
        

# mmpdb drop_index <filename>
def drop_index_command(parser, args):
    mmpdb = dbutils.open_database_from_args_or_exit(args)
    with mmpdb.atomic():
        schema.drop_index(mmpdb)
    mmpdb.execute("VACUUM")

##### loadprops/reaggregate

def reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                           cursor, reporter):
    # Mapping from rule environment id to rule environment statistics id
    
    reporter.update("Computing aggregate statistics")
    num_pairs = dataset.get_num_pairs(cursor=cursor)
    
    all_pairs = dataset.iter_pairs(cursor=cursor)
    all_pairs = reporter.progress(all_pairs, "Computing aggregate statistics", num_pairs)

    stats_info = []
    seen_rule_environment_ids = set()
    for rule_environment_id, rule_environment_pairs in itertools.groupby(
            all_pairs, (lambda pair: pair.rule_environment_id)):

        seen_rule_environment_ids.add(rule_environment_id)
        rule_environment_pairs = list(rule_environment_pairs)  # now a list, not iterator
        
        for property_name_id in property_name_ids:
            deltas = []
            compound_values = compound_values_for_property_name_id[property_name_id]
            for pair in rule_environment_pairs:
                value1 = compound_values.get(pair.compound1_id, None)
                if value1 is None:
                    continue
                value2 = compound_values.get(pair.compound2_id, None)
                if value2 is None:
                    continue
                deltas.append(value2-value1)
            if deltas:
                stats = index_algorithm.compute_aggregate_values(deltas)
                stats_info.append( (rule_environment_id, property_name_id, stats) )

    # Need to figure out if the statistics exist or need to be created
    reporter.report("Generated %d rule statistics (%d rule environments, %d properties)"
                    % (len(stats_info), len(seen_rule_environment_ids), len(property_name_ids)))
    reporter.update("Getting information about which rule statistics exist...")
    existing_stats_ids = dataset.get_rule_environment_statistics_mapping(
        property_name_ids, cursor=cursor)

    stats_info_progress = reporter.progress(
        stats_info, "Updating statistics table", len(stats_info))
    seen_stats_ids = set()
    num_updated = num_added = 0
    for (rule_environment_id, property_name_id, stats) in stats_info_progress:
        key = (rule_environment_id, property_name_id)
        stats_id = existing_stats_ids.get(key, None)
        if stats_id is not None:
            dataset.update_rule_environment_statistics(stats_id, stats)
            seen_stats_ids.add(stats_id)
            num_updated += 1
        else:
            dataset.add_rule_environment_statistics(rule_environment_id, property_name_id, stats)
            num_added += 1

    to_delete = set(existing_stats_ids.values()) - seen_stats_ids
    num_deleted = len(to_delete)
    if to_delete:
        delete_report = reporter.progress(
            iter(to_delete), "Deleting rule statistics", num_deleted)
        while 1:
            ids = list(itertools.islice(delete_report, 0, 1000))
            if not ids:
                break
            dataset.delete_rule_environment_statistics(ids)

    reporter.report("Number of rule statistics added: %d updated: %d deleted: %d"
                    % (num_added, num_updated, num_deleted))
        
        
# mmpdb loadprops <filename> -p <properties_filename>
def loadprops_command(parser, args):
    from . import properties_io
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()
    reporter.report("Using dataset: %s" % (dataset.title,))
        
    if args.properties is None:
        reporter.report("Reading properties from stdin")
        properties_file = sys.stdin
        close = None
        source = "<stdin>"
    else:
        reporter.report("Reading properties from %r" % (args.properties,))
        try:
            properties_file = open(args.properties, "U")
        except IOError as err:
            parser.error("Cannot open properties file: %s" % (err,))
        close = properties_file.close
        source = args.properties

    try:
        try:
            with properties_file:
                properties = properties_io.load_properties(properties_file, reporter)
        except ValueError as err:
            parser.error("Problem reading --properties file %r: %s"
                         % (args.properties, err))
    finally:
        if close is not None:
            close()

    reporter.report("Read %d properties for %d compounds from %r"
                    % (len(properties.property_columns), len(properties.id_column),
                       source))

    public_id_to_id = dataset.get_public_id_to_id_table(c)
    
    compound_ids = [public_id_to_id.get(public_id, None) for public_id in properties.id_column]
    num_missing = compound_ids.count(None)
    if num_missing:
        reporter.report("%d compounds from %r are not in the dataset at %r"
                        % (num_missing, source, args.databases[0]))
        ## missing = [public_id for public_id in properties.id_column if public_id not in public_id_to_id]
        ## del missing[6:]
        ## if len(missing) > 5:
        ##     reporter.warning("First five missing records: %r" % (missing[:5],))
        ## else:
        ##     reporter.warning("Missing records: %r" % (missing,))
            

    UPDATE_COMPOUND_PROPERTY_SQL = db.SQL(
        "UPDATE compound_property "
        "      SET value = ?"
        "       WHERE compound_id = ?"
        "         AND property_name_id = ?")
    INSERT_COMPOUND_PROPERTY_SQL = db.SQL(
        "INSERT INTO compound_property (compound_id, property_name_id, value) "
        " VALUES (?, ?, ?)")
    
    with db.atomic():
        # Remember which compound properties exist, so I can tell if a
        # value should replace an existing value or create a new value.
        c.execute(db.SQL(
            "SELECT compound_id, property_name_id from compound_property"))
        seen_properties = dict((key, False) for key in c)

        compound_values_for_property_name_id = {}
        property_name_ids = []
        
        for property_name, property_values in properties.iter_properties():
            property_name_id = dataset.get_or_add_property_name(property_name)
            property_name_ids.append(property_name_id)
            #reporter.report("Loading property %r (id %d)" % (property_name.name, property_name.id))

            num_created = num_updated = 0
            compound_values_for_property_name_id[property_name_id] = compound_values = {}
            
            for compound_id, value in zip(compound_ids, property_values):
                if compound_id is None:
                    # Property specified but not in database
                    continue
                if value is None:
                    # Property value is "*", meaning it isn't available
                    continue
                key = (compound_id, property_name_id)
                if key in seen_properties:
                    seen_properties[key] = True
                    num_updated += 1
                    c.execute(UPDATE_COMPOUND_PROPERTY_SQL,
                              (value, compound_id, property_name_id))
                else:
                    num_created += 1
                    c.execute(INSERT_COMPOUND_PROPERTY_SQL,
                              (compound_id, property_name_id, value))
                compound_values[compound_id] = value
            reporter.report(
                "Imported %d %r records (%d new, %d updated)."
                % (num_updated + num_created, property_name, num_created, num_updated))

        # Remove existing compound properties where the property name was in the
        # properties file but the where the file did not specify a value.
        properties_to_delete = [key for key, was_updated in seen_properties.items()
                                if not was_updated and key[1] in property_name_ids]
        if properties_to_delete:
            dataset.delete_compound_properties(properties_to_delete)
        
        reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                               cursor=c, reporter=reporter)

        # Check if any of the properties are completely gone
        if properties_to_delete:
            for property_name_id in property_name_ids:
                n = dataset.get_num_compound_properties(property_name_id, cursor=c)
                if n == 0:
                    dataset.delete_property_name_id(property_name_id, cursor=c)

        # Update the environment statistics
        reporter.update("Updating environment statistics count ...")
        num = schema._get_one(c.execute("SELECT count(*) from rule_environment_statistics"))
        c.execute("UPDATE dataset SET num_rule_environment_stats=?", (num,))
        
        reporter.update("Commiting changed ...")
                
    reporter.report("Loaded all properties and re-computed all rule statistics.")

            
def reaggregate_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()
    
    with db.atomic():
        # Which properties to reaggregate?
        property_names = command_support.get_property_names_or_error(parser, args, dataset)
        reporter.report("Reaggregating %d properties" % (len(property_names),))

        # Convert them into ids
        property_name_ids = [dataset.get_property_name_id(name) for name in property_names]

        # Get their values
        compound_values_for_property_name_id = dict(
            (property_name_id, dataset.get_property_values(property_name_id))
             for property_name_id in property_name_ids)

        # Reaggregate
        reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                               cursor=c, reporter=reporter)
                

# mmpdb smicat [--input-smiles] [-o filename] filename
def smicat_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()
    input_smiles = args.input_smiles

    output_filename = args.output
    with fileio.open_output(output_filename, output_filename) as outfile:
        for compound_row in dataset.iter_compounds(cursor=c):
            if input_smiles:
                smiles = compound_row.input_smiles
            else:
                smiles = compound_row.clean_smiles
            outfile.write("%s %s\n" % (smiles, compound_row.public_id))

# mmpdb propcat [--property <name>]* [--all] [-o filename] filename
def propcat_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()

    show_all = args.all

    property_names = command_support.get_property_names_or_error(parser, args, dataset)

    property_values_list = []
    for property_name in property_names:
        property_name_id = dataset.get_property_name_id(property_name)
        property_values_list.append(dataset.get_property_values(property_name_id))

    output_filename = args.output
    with fileio.open_output(output_filename, output_filename) as outfile:
        print("ID", *property_names, sep="\t", file=outfile)
        for compound_row in dataset.iter_compounds(cursor=c):
            columns = [compound_row.public_id]
            is_empty = True
            for property_values in property_values_list:
                value = property_values.get(compound_row.id, None)
                if value is None:
                    columns.append("*")
                else:
                    columns.append(value)
                    is_empty = False
            if show_all or not is_empty:
                print(*columns, sep="\t", file=outfile)
        


# mmpdb rulecat? or would paircat be more useful? both?
