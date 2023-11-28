"""Helper functions related to fragmentation"""

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

from .click_utils import (
    IntChoice,
    die,
    pop_known_args,
    positive_int_or_none,
    nonnegative_int,
    set_click_attrs,
)

from .. import config
from .. import smarts_aliases


def callback_chain(*funcs):
    def caller(ctx, param, value):
        for func in funcs:
            value = func(ctx, param, value)
        return value

    return caller


EXCLUSION_CUTS = "mmpdblib.exclusion.cuts"


def mutual_exclusion_cuts(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return value
    name = param.name
    prev_name = ctx.meta.get(EXCLUSION_CUTS, None)
    if prev_name is None:
        ctx.meta[EXCLUSION_CUTS] = name
        return value
    if prev_name == name:
        return value
    raise click.UsageError("Cannot specify more than one of --cut-smarts, --cut-rgroup, or --cut-rgroup-file")


# Record information about which fragmentation options are used so
# "fragment" can detect conflicts in mixing these with --cache
FRAGMENTATION_OPTION_NAMES = "mmpdblib.fragment_option_names"


def record_used_options(ctx, param, value):
    used_names = ctx.meta.setdefault(FRAGMENTATION_OPTION_NAMES, set())
    if not value:
        return value
    used_names.add(param.name)
    return value


#### @add_fragment_options


def add_fragment_options(command):
    OPTS = config.DEFAULT_FRAGMENT_OPTIONS

    param_names = []

    def add_option(*args, **kwargs):
        # Keep track of the parameter names used
        param_names.append(args[0].lstrip("-").replace("-", "_"))

        # Want to track which fragment options were used, in case one is
        # specified with --cache.  This only works if the default is None,
        # so we can't use
        #    default = OPT.option_name
        # These defaults are instead set at the very end.
        callback = record_used_options
        if "callback" in kwargs:
            callback = callback_chain(callback, kwargs["callback"])
        kwargs["callback"] = callback

        click.option(*args, **kwargs)(command)

    add_option(
        "--max-heavies",
        type=positive_int_or_none(),
        help=f"Maximum number of non-hydrogen atoms, or 'none' (default: {OPTS.max_heavies})",
    )

    add_option(
        "--max-rotatable-bonds",
        type=positive_int_or_none(),
        help=f"Maximum number of rotatable bonds (default: {OPTS.max_rotatable_bonds})",
    )

    add_option(
        "--rotatable-smarts",
        metavar="SMARTS",
        help=f"SMARTS pattern to detect rotatable bonds (default: {OPTS.rotatable_smarts!r})",
    )

    add_option(
        "--salt-remover",
        metavar="FILENAME",
        help=(
            "File containing RDKit SaltRemover definitions. The default ('<default>') "
            "uses RDKit's standard salt remover. Use '<none>' to not remove salts."
        ),
    )

    ## These are mutually exclusive group
    alias_names = ", ".join(repr(alias.name) for alias in smarts_aliases.cut_smarts_aliases)

    add_option(
        "--cut-smarts",
        metavar="SMARTS",
        help=(
            f"Alternate SMARTS pattern to use for cutting (default: {OPTS.cut_smarts!r}), "
            f"or use one of: {alias_names}"
        ),
        callback=mutual_exclusion_cuts,
    )

    add_option(
        "--cut-rgroup",
        metavar="SMILES",
        multiple=True,
        help="Cut on the attachment point for the given R-group SMILES",
        callback=mutual_exclusion_cuts,
    )

    add_option(
        "--cut-rgroup-file",
        metavar="FILENAME",
        help="Read R-group SMILES from the named file",
        callback=mutual_exclusion_cuts,
    )

    ##
    add_option(
        "--num-cuts",
        type=IntChoice(["1", "2", "3"]),
        help=f"Number of cuts to use (default: {OPTS.num_cuts})",
    )

    add_option(
        "--min-heavies-per-const-frag",
        type=nonnegative_int(),
        metavar="N",
        help=(
            "Ignore fragmentations where one or more constant fragments have "
            f"fewer than N heavy atoms (default: {OPTS.min_heavies_per_const_frag})"
        ),
    )
    
    add_option(
        "--min-heavies-total-const-frag",
        type=nonnegative_int(),
        metavar="N",
        help=(
            "Ignore fragmentations where there are fewer than N heavy atoms in the "
            "total constant fragment  (default: {OPTS.min_heavies_total_const_frag})"
        ),
    )
    
    add_option(
        "--max-up-enumerations",
        type=nonnegative_int(),
        metavar="N",
        help=(
            "Maximum number of up-enumerations "
            f"(default: {OPTS.max_up_enumerations})"
        ),
    )

    # Wrap the command to convert the fragment option parameters
    # into a single object
    def make_fragment_options_wrapper(**kwargs):
        # Fill in the defaults, or use None if there aren't defaults (eg, for
        # --cut-rgroups and --cut-rgroup-file).
        popped_kwargs = pop_known_args(param_names, kwargs, OPTS)

        kwargs["fragment_options"] = make_fragment_options(**popped_kwargs)

        # Forward to the command
        return command(**kwargs)

    set_click_attrs(make_fragment_options_wrapper, command)

    return make_fragment_options_wrapper


########


def make_fragment_options(
    *,
    max_heavies,
    max_rotatable_bonds,
    rotatable_smarts,
    cut_rgroup_file,
    cut_rgroup,
    cut_smarts,
    num_cuts,
    salt_remover,
    min_heavies_per_const_frag,
    min_heavies_total_const_frag,
    max_up_enumerations,
):
    from .. import (
        fragment_types,
        rgroup2smarts,
    )

    if cut_rgroup_file is not None:
        try:
            cut_smarts = rgroup2smarts.get_recursive_smarts_from_cut_filename(cut_rgroup_file)
        except OSError as err:
            die(f"Cannot use --cut-rgroup-file: {cut_rgroup_file!r}: {err}")

        except rgroup2smarts.ParseError as err:
            die(f"Cannot parse --cut-rgroup-file: {cut_rgroup_file!r}: {err}")

        except rgroup2smarts.ConversionError as err:
            die(f"Error in --cut-rgroup-file: {cut_rgroup_file!r}: {err}")

    elif cut_rgroup:
        try:
            cut_smarts = rgroup2smarts.get_recursive_smarts_from_cut_rgroups(
                cut_rgroup, source="--cut-rgroup", offset=1
            )
        except rgroup2smarts.ConversionError as err:
            die(str(err))

    else:
        # Resolve any alias
        if cut_smarts in smarts_aliases.cut_smarts_aliases_by_name:
            cut_smarts = smarts_aliases.cut_smarts_aliases_by_name[cut_smarts].smarts

    method = "chiral"

    if max_heavies == "none":
        max_heavies = None

    if max_rotatable_bonds == "none":
        max_rotatable_bonds = None

    return fragment_types.FragmentOptions(
        max_heavies=max_heavies,
        max_rotatable_bonds=max_rotatable_bonds,
        rotatable_smarts=rotatable_smarts,
        cut_smarts=cut_smarts,
        num_cuts=num_cuts,
        salt_remover=salt_remover,
        method=method,
        min_heavies_per_const_frag=min_heavies_per_const_frag,
        min_heavies_total_const_frag=min_heavies_total_const_frag,
        max_up_enumerations=max_up_enumerations,
    )
