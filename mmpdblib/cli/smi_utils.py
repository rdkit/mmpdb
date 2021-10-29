"Utilities for working with SMILES file input reading"

import click

class SmiInputOptions:
    def __init__(self, format, delimiter, has_header):
        self.format = format
        self.delimiter = delimiter
        self.has_header = has_header
        

def add_input_options(f):
    def add_option(*args, **kwargs):
        click.option(*args, **kwargs)(f)

    add_option(
        "--in",
        "-i",
        "in_format",
        type = click.Choice(["smi", "smi.gz"]),
        help = (
            "input structuture format (one of 'smi', 'smi.gz'). "
            "If not specified, use the filename extension or default to 'smi'."
            ),
        )
    
    add_option(
        "--delimiter",
        default="whitespace",
        type = click.Choice(["whitespace", "to-eol", "comma", "tab", "space", "native"]),
        help = (
            "SMILES file delimiter style (one of 'whitespace' (default), 'to-eol', "
            "'comma', 'tab', 'space', or 'native')"
            ),
        )
    
    add_option(
        "--has-header",
        default = False,
        help = "skip the first line, which is the header line",
        )

    # Wrap the command to convert the fragment option parameters
    # into a single object
    
    def make_input_options_wrapper(**kwargs):
        kwargs["input_options"] = SmiInputOptions(
            format = kwargs.pop("in_format"),
            
            delimiter = kwargs.pop("delimiter"),
            has_header = kwargs.pop("has_header")
            )
        return f(**kwargs)

    make_input_options_wrapper.__name__ = f.__name__
    make_input_options_wrapper.__click_params__ = f.__click_params__
            
    return make_input_options_wrapper
    
