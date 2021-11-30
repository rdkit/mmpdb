"Implement the 'loadprops' command"

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

import sys
import click

from .click_utils import (
    command,
    die,
    open_database_from_options_or_exit,
    add_single_database_parameters,
)


loadprops_epilog = """
Load structure property values from a CSV into a data set.

The property file contains a data table. Here is an example:

\b
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

\b
  % mmpdb loadprops --properties MP.csv mmpdb.db
  
"""


@command(epilog=loadprops_epilog)
@click.option(
    "--properties",
    "-p",
    "properties_filename",
    metavar="FILENAME",
    help="File containing the identifiers to use and optional physical properties",
)
@add_single_database_parameters()
@click.pass_obj
def loadprops(
    reporter,
    properties_filename,
    database_options,
):
    """Load properties for existing structures

    DATABASE: the mmpdb database to update with the new properties
    """
    from .. import properties_io, dbutils, schema

    db = open_database_from_options_or_exit(database_options)
    c = db.get_cursor()
    dataset = db.get_dataset()
    reporter.report(f"Using dataset: {dataset.title}")

    if properties_filename is None:
        reporter.report("Reading properties from stdin")
        properties_file = sys.stdin
        close = None
        source = "<stdin>"
    else:
        reporter.report(f"Reading properties from {properties_filename!r}")
        try:
            properties_file = open(properties_filename)
        except IOError as err:
            die(f"Cannot open properties file: {err}")
        close = properties_file.close
        source = properties_filename

    try:
        try:
            with properties_file:
                properties = properties_io.load_properties(properties_file, reporter)
        except ValueError as err:
            die(f"Problem reading --properties file {properties_filename}: {err}")
    finally:
        if close is not None:
            close()

    reporter.report(
        f"Read {len(properties.property_columns)} properties for "
        f"{len(properties.id_column)} compounds "
        f"from {source!r}"
    )
    public_id_to_id = dataset.get_public_id_to_id_table(c)

    compound_ids = [public_id_to_id.get(public_id, None) for public_id in properties.id_column]
    num_missing = compound_ids.count(None)
    if num_missing:
        reporter.report(
            f"{num_missing} compounds from {source!r} are not in the dataset at {database_options.database!r}"
        )
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
        "         AND property_name_id = ?"
    )
    INSERT_COMPOUND_PROPERTY_SQL = db.SQL(
        "INSERT INTO compound_property (compound_id, property_name_id, value) " " VALUES (?, ?, ?)"
    )

    with db.atomic():
        # Remember which compound properties exist, so I can tell if a
        # value should replace an existing value or create a new value.
        c.execute(db.SQL("SELECT compound_id, property_name_id from compound_property"))
        seen_properties = dict((key, False) for key in c)

        compound_values_for_property_name_id = {}
        property_name_ids = []

        for property_name, property_values in properties.iter_properties():
            property_name_id = dataset.get_or_add_property_name(property_name)
            property_name_ids.append(property_name_id)
            # reporter.report("Loading property %r (id %d)" % (property_name.name, property_name.id))

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
                    c.execute(
                        UPDATE_COMPOUND_PROPERTY_SQL,
                        (value, compound_id, property_name_id),
                    )
                else:
                    num_created += 1
                    c.execute(
                        INSERT_COMPOUND_PROPERTY_SQL,
                        (compound_id, property_name_id, value),
                    )
                compound_values[compound_id] = value
            reporter.report(
                f"Imported {num_updated + num_created} {property_name!r} records "
                f"({num_created} new, {num_updated} updated)."
            )

        # Remove existing compound properties where the property name was in the
        # properties file but the where the file did not specify a value.
        properties_to_delete = [
            key for key, was_updated in seen_properties.items() if not was_updated and key[1] in property_name_ids
        ]
        if properties_to_delete:
            dataset.delete_compound_properties(properties_to_delete)

        dbutils.reaggregate_properties(
            dataset,
            property_name_ids,
            compound_values_for_property_name_id,
            cursor=c,
            reporter=reporter,
        )

        # Check if any of the properties are completely gone
        if properties_to_delete:
            for property_name_id in property_name_ids:
                n = dataset.get_num_compound_properties(property_name_id, cursor=c)
                if n == 0:
                    dataset.delete_property_name_id(property_name_id, cursor=c)

        # Update the environment statistics
        reporter.update("Updating environment statistics count ...")
        c.execute("SELECT count(*) from rule_environment_statistics")
        num = schema._get_one(c)
        c.execute(db.SQL("UPDATE dataset SET num_rule_environment_stats=?"), (num,))

        reporter.update("Commiting changed ...")

    reporter.report("Loaded all properties and re-computed all rule statistics.")
