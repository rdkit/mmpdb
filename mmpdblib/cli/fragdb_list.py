import sys
import dataclasses
import click

from .click_utils import (
    command,
    add_multiple_databases_parameters,
    )

@dataclasses.dataclass
class FragDBInfo:
    filename: str
    num_compounds: int
    num_error_compounds: int
    num_fragmentations: int
    num_constants: int
    num_variables: int
    max_num_pairs: int
    options: object

    def get_cols(self):
        return [
            self.filename,
            str(self.num_compounds),
            str(self.num_error_compounds),
            str(self.num_fragmentations),
            str(self.num_constants),
            str(self.num_variables),
            str(self.max_num_pairs),
            ]

def write(terms):
    sys.stdout.write(" ".join(terms) + "\n")

def get_info(filename, reporter):
    from .. import fragment_db
    try:
        db = fragment_db.open_fragdb(filename)
    except IOError as err:
        reporter.warning(f"Cannot open database: {err} -- Skipping.")
        return None
    except ValueError as err:
        reporter.warning(f"Cannot use database: {err} -- Skipping.")
        return None

    def _get_one(sql):
        c.execute(sql)
        for (n,) in c:
            return n
        raise AssertionError("cannot get one", sql)
    
    with db:
        c = db.cursor()
        return FragDBInfo(
            filename = filename,
            num_compounds = _get_one(
                "SELECT COUNT(*) FROM record"
                ),
            num_error_compounds = _get_one(
                "SELECT COUNT(*) FROM error_record",
                ),
            num_fragmentations = _get_one(
                "SELECT COUNT(*) FROM fragmentation",
                ),
            num_constants = _get_one(
                "SELECT COUNT(DISTINCT constant_smiles) FROM fragmentation",
                ),
            num_variables = _get_one(
                "SELECT COUNT(DISTINCT variable_smiles) FROM fragmentation",
                ),
            max_num_pairs = _get_one(
                # XXX Should I have a +1 for single-cut constants?
                "SELECT SUM(i*(i-1)/2) FROM (SELECT COUNT(*) AS i FROM fragmentation GROUP BY constant_smiles)",
                ),
            options = db.options,
            )
            
    
@command(
    name = "fragdb_list",
    )

@click.option(
    "--all",
    "-a",
    "show_all",
    is_flag = True,
    default = False,
    help = "Include option information",
    )

@add_multiple_databases_parameters()
@click.pass_obj
def fragdb_list(
        reporter,
        databases_options,
        show_all,
        ):
    """Summarize zero or more fragdb databases

    If no DATABASE is given then look for '*.fragdb' in the current directory.
    """

    databases = databases_options.databases
    if not databases:
        import glob
        databases = glob.glob("*.fragdb")
        databases.sort()

    col_headers = [        
        "Name",
        "#recs",
        "#errs",
        "#frags",
        "#consts",
        "#vars",
        "max.#pairs",
        ]
    col_sizes = [len(s) for s in col_headers]
        
    info_list = []
    rows = []

    for database in databases:
        info = get_info(database, reporter)
        if info is None:
            continue
        info_list.append(info)
        
        cols = info.get_cols()
        rows.append(cols)

        # Figure out the columns sizes
        for i, col in enumerate(cols):
            n = len(col)
            if n > col_sizes[i]:
                col_sizes[i] = n

    write(header.center(col_size)
              for header, col_size in zip(col_headers, col_sizes))

    for info, row in zip(info_list, rows):
        write(col.rjust(col_size)
                  for col, col_size in zip(row, col_sizes))
        if show_all:
            sys.stdout.write("        Fragment options:\n")
            d = info.options.to_dict()
            for k, v in d.items():
                sys.stdout.write(f"          {k}: {v}\n")
