"""Helper functions related to using click"""

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
    def __init__(self):
        super().__init__(0)
    
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
        click.echo(msg, error=True)
    raise SystemExit(1)
