"""Helper functions related to using click"""

# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2021, Andrew Dalke Scientific AB
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

from .. import config
from .. import smarts_aliases

class FormatEpilog:
    def format_epilog(self, ctx, formatter):
        if self.epilog:
            formatter.write_paragraph()
            # Change the indentation level so it isn't
            # quite at the "Commands:" level.
            old_value = formatter.indent_increment
            formatter.indent_increment = 1
            try:
                formatter.indent()
                formatter.write_text(self.epilog)
            finally:
                formatter.dedent()
                formatter.indent_increment = old_value

# Have the command order match the insertion order.
# https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
class OrderedGroup(FormatEpilog, click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(OrderedGroup, self).__init__(name, commands, **attrs)
        # the registered subcommands by their exported names.
        self.commands = commands or {}

    def list_commands(self, ctx):
        return self.commands

class Command(FormatEpilog, click.Command):
    pass

def command(**kwargs):
    d = {"cls": Command}
    d.update(kwargs)
    return click.command(**d)

def ordered_group(**kwargs):
    d = {"cls": OrderedGroup}
    d.update(kwargs)
    return click.group(**d)

#### ParamType


class positive_int_or_none(click.ParamType):
    name = "N"
    def convert(self, value, param, ctx):
        if value is None or value == "none":
            return None

        msg = "must be a positive integer or 'none'"
        if isinstance(value, str):
            try:
                value = int(value)
            except ValueError:
                self.fail(msg)

        if not isinstance(value, int):
            self.fail(msg)

        if value <= 0:
            self.fail(msg)
            
        return value


class IntChoice(click.Choice):
    def convert(self, value, param, ctx):
        if isinstance(value, int):
            value = str(value)
        return int(super().convert(value, param, ctx))
        
class nonnegative_float(click.ParamType):
    name = "FLT"
    def convert(self, value, param, ctx):
        msg = "must be a positive float or zero"
        if isinstance(value, str):
            try:
                value = float(value)
            except ValueError:
                self.fail(msg)
        if not (value >= 0.0):
            self.fail(msg)
        return value
            
class positive_float(click.ParamType):
    name = "FLT"
    def convert(self, value, param, ctx):
        msg = "must be a positive float"
        if isinstance(value, str):
            try:
                value = float(value)
            except ValueError:
                self.fail(msg)
        if not (value > 0.0):
            self.fail(msg)
        return value
            
class nonnegative_int(click.IntRange):
    name = "N"
    def __init__(self):
        super().__init__(0)

class positive_int(click.IntRange):
    name = "N"
    def __init__(self):
        super().__init__(1)

class radius_type(IntChoice):
    name = "R"
    def __init__(self):
        super().__init__(["0", "1", "2", "3", "4", "5"])
    
def name_to_command_line(s):
    return "--" + s.replace("_", "-")
    
def set_click_attrs(dest, src):
    dest.__name__ = src.__name__
    dest.__doc__ = src.__doc__
    dest.__click_params__ = src.__click_params__

# Used in the click call wrappers
def pop_known_args(names, kwargs, opts):
    popped_kwargs = {}
    for name in names:
        value = kwargs.pop(name)
        if value is None:
            value = getattr(opts, name, None)

        popped_kwargs[name] = value
    return popped_kwargs
    
def die(*msgs):
    for msg in msgs:
        click.echo(msg, err=True)
    raise SystemExit(1)


###### Shared options

class DatabaseOptions:
    def __init__(self, database, copy_to_memory = False):
        self.database = database
        self.copy_to_memory = copy_to_memory

class DatabasesOptions:
    def __init__(self, databases, copy_to_memory = False):
        self.databases = databases
        self.copy_to_memory = copy_to_memory

        
def _in_memory_option(command):
    click.option(
        "--in-memory",
        is_flag=True,
        default=False,
        help = "load the SQLite database into memory before use (requires APSW)"
        )(command)
        
def add_single_database_parameters(add_in_memory = False):
    def add_single_database_parameters_decorator(command):
        click.argument(
            "database",
            metavar = "DATABASE",
            )(command)
        
        if add_in_memory:
            _in_memory_option(command)
            
        def wrapped_command(**kwargs):
            popped_kwargs = {
                "database":  kwargs.pop("database")
                }
            if add_in_memory:
                popped_kwargs["copy_to_memory"] = kwargs.pop("in_memory")
            kwargs["database_options"] = DatabaseOptions(**popped_kwargs)
            return command(**kwargs)

        set_click_attrs(wrapped_command, command)
        return wrapped_command
    return add_single_database_parameters_decorator

def add_multiple_databases_parameters(add_in_memory = False):
    def add_multiple_databases_parameters_decorator(command):
        click.argument(
            "database",
            nargs = -1,
            metavar="DATABASE",
            )(command)
        
        if add_in_memory:
            _in_memory_option(command)
            
        def wrapped_command(**kwargs):
            popped_kwargs = {
                "databases": kwargs.pop("database")  # is 0 or more
                }
            if add_in_memory:
                popped_kwargs["copy_to_memory"] = kwargs.pop("in_memory")
            kwargs["databases_options"] = DatabasesOptions(**popped_kwargs)
            return command(**kwargs)

        set_click_attrs(wrapped_command, command)
        return wrapped_command
    
    return add_multiple_databases_parameters_decorator

def add_single_property(command):
    click.option(
        "--property",
        "-p",
        "property_name",
        metavar="NAME",
        required = True,
        help = "property to use",
        )(command)
    return command

def add_multiple_properties(command):
    click.option(
        "--property",
        "-p",
        "property_names",
        metavar = "NAME",
        multiple = True,
        help = "property to use (may be specified multiple times)",
        )(command)
    click.option(
        "--no-properties",
        is_flag = True,
        default = False,
        help = "don't use any properties",
        )(command)
    return command

## Rule selection

class parse_where(click.ParamType):
    name = "EXPR"
    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value
        from .. import analysis_algorithms
        try:
            return analysis_algorithms.get_where_function(value)
            
        except ValueError as err:
            self.fail(str(err), param, ctx)
        
        except analysis_algorithms.EvalError as err:
            self.fail(str(err), param, ctx)

class parse_score(click.ParamType):
    name = "EXPR"
    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            return value
        from .. import analysis_algorithms
        try:
            return analysis_algorithms.get_score_function(value)
        except ValueError as err:
            self.fail(str(err), param, ctx)
        
        except analysis_algorithms.EvalError as err:
            self.fail(str(err), param, ctx)



class parse_cutoff_list(click.ParamType):
    name = "LIST"
    def convert(self, value_obj, param, ctx):
        if isinstance(value_obj, (list, tuple)):
            return value_obj
        
        prev = None
        values = []
        for term in value_obj.split(","):
            try:
                value = int(term)
            except ValueError as err:
                self.fail(
                    f"could not parse {term} as an integer: {err}",
                    param,
                    ctx,
                    )

            if value < 0:
                self.fail(
                    "threshold values must be non-negative",
                    param,
                    ctx,
                    )

            if prev is not None and prev <= value:
                self.fail(
                    "threshold values must be in decreasing order",
                    param,
                    ctx,
                    )
                    
            
            prev = value
            values.append(value)

        if not values:  # Let people specify ""
            return [0]

        return values        

def add_rule_selection_options(command):
    OPTS = config.DEFAULT_RULE_SELECTION_OPTIONS

    def add_option(*args, **kwargs):
        click.option(*args, **kwargs)(command)

    add_option(
        "--where",
        default=OPTS.where,
        type = parse_where(),
        help = "select only rules for which the expression is true",
        )

    add_option(
        "--score",
        default= OPTS.score,
        type = parse_score(),
        help="use to break ties when multiple rules produce the same SMILES",
        )

    msg = ",".join(str(i) for i in OPTS.cutoff_list)
    add_option(
        "--rule-selection-cutoffs",
        default = OPTS.cutoff_list,
        type = parse_cutoff_list(),
        help=(
            "evaluate rule environments with the given minimum pair count. If multiple "
            "counts are given, consider them in turn until there is a selected environment. "
            f"(default: '{msg}')"
            ),
        )

    def wrapped_command(**kwargs):
        kwargs["rule_selection_options"] = config.RuleSelectionOptions(
            where = kwargs.pop("where"),
            score = kwargs.pop("score"),
            cutoff_list = kwargs.pop("rule_selection_cutoffs"),
            )
        return command(**kwargs)

    set_click_attrs(wrapped_command, command)
    return wrapped_command

## Properties

def get_property_names_or_error(dataset, *, property_names, no_properties=False, all_properties=False):
    property_names = []

    if property_names:
        if no_properties:
            raise click.UsageError("Cannot specify --property and --no-properties")
        if all_properties:
            raise click.UsageError("Cannot specify --property and --all-properties")
    
    known_names = dataset.get_property_names()
    if not property_names:
        if no_properties:
            return []
        # Use all of the properties
        return known_names

    known_names = set(known_names)
    seen = set()
    for name in property_names:
        # If a property is specified multiple times, only use the first one.
        if name in seen:
            continue
        seen.add(name)

        if name not in known_names:
            die(f"--property {name!r} is not present in the database")
        property_names.append(name)
    return property_names
