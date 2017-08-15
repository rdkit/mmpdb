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

import sys
import multiprocessing
import signal

from rdkit import Chem
from rdkit.Chem import SaltRemover


from . import config
from . import command_support
from . import reporters
from . import fileio
from . import fragment_algorithm
from . import fragment_io
from . import fragment_types
from . import smarts_aliases

from .fileio import FileFormatError

###


def parse_rotatable_smarts(smarts):
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError("unable to parse SMARTS")
    if pattern.GetNumAtoms() != 2:
        raise ValueError("rotatable SMARTS must match exactly two atoms")
    if pattern.GetNumBonds() != 1:
        raise ValueError("rotatable SMARTS must connect both atoms")
    
    return pattern

def parse_cut_smarts(smarts):
    if smarts in smarts_aliases.cut_smarts_aliases_by_name:
        smarts = smarts_aliases.cut_smarts_aliases_by_name[smarts].smarts
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError("unable to parse SMARTS")
    if pattern.GetNumAtoms() != 2:
        raise ValueError("cut SMARTS must match exactly two atoms")
    if pattern.GetNumBonds() != 1:
        raise ValueError("cut SMARTS must connect both atoms")
    return pattern


def parse_salt_remover(salt_remover_filename):
    if salt_remover_filename == "<none>":
        # Don't use a salt remover
        return None
    elif salt_remover_filename == "<default>":
        remover = SaltRemover.SaltRemover()
    else:
        # raises an I/O error on failure
        remover = SaltRemover.SaltRemover(salt_remover_filename)

    # Work-around rdkit issue #541, where SMARTS failures get
    # turned into None patterns instead of raising an exception
    # during SaltRemover initialization.
    for i, pattern in enumerate(remover.salts, 1):
        if pattern is None:
            raise ValueError("SMARTS pattern #%d is invalid" % (i,))
    return remover

def parse_num_cuts(num_cuts):
    if num_cuts in (1, 2, 3):
        return num_cuts
    raise ValueError("must be 1, 2, or 3")

def parse_method(method):
    if method in ("", None, "chiral"):
        return fragment_algorithm.fragment_mol
    raise ValueError("must be 'chiral'")


###

def get_fragment_options_from_args(parser, args):
    specified_args = set()
    def get(name):
        value = getattr(args, name)
        if value is None:
            return getattr(config.DEFAULT_FRAGMENT_OPTIONS, name)
        specified_args.add(name)
        return value

    max_heavies = get("max_heavies")
    if max_heavies == "none":
        max_heavies = None
        
    max_rotatable_bonds = get("max_rotatable_bonds")
    if max_rotatable_bonds == "none":
        max_rotatable_bonds = None

    rotatable_smarts = get("rotatable_smarts")

    cut_smarts = get("cut_smarts")
    # Resolve any alias
    if cut_smarts in smarts_aliases.cut_smarts_aliases_by_name:
        cut_smarts = smarts_aliases.cut_smarts_aliases_by_name[cut_smarts].smarts
    
    num_cuts = get("num_cuts")
    assert num_cuts in (1, 2, 3)

    salt_remover = get("salt_remover")

    #method = get("method")
    method = "chiral"

    return specified_args, config.FragmentOptions(
        max_heavies = max_heavies,
        max_rotatable_bonds = max_rotatable_bonds,
        rotatable_smarts =rotatable_smarts,
        cut_smarts = cut_smarts,
        num_cuts = num_cuts,
        salt_remover = salt_remover,
        method = method,
        )
###

class FragmentFilter(object):
    def __init__(self, max_heavies, max_rotatable_bonds,
                 rotatable_pattern, salt_remover,
                 cut_pattern, num_cuts, method,
                 options
                 ):
        if rotatable_pattern is None:
            max_rotatable_bonds = None
            
        self.max_heavies = max_heavies
        self.max_rotatable_bonds = max_rotatable_bonds
        self.rotatable_pattern = rotatable_pattern
        self.salt_remover = salt_remover
        self.cut_pattern = cut_pattern
        self.num_cuts = num_cuts
        self.method = method
        self.options = options

    def normalize(self, mol):
        # XXX Remove existing isotope labels?
        if self.salt_remover is not None:
            desalted_mol = self.salt_remover.StripMol(mol)
            if not mol.GetNumAtoms():
                return ("no non-salts", desalted_mol)
        else:
            desalted_mol = mol

        return None, desalted_mol
        
    def apply_filters(self, mol):
        num_heavies = 0
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num <= 1:
                if atomic_num == 0:
                    return "contains '*' atoms"
                else:
                    # should I worry about explicit hydrogens?
                    pass
            else:
                num_heavies += 1

        if num_heavies < 2:
            return "not enough heavy atoms"
                
        max_heavies = self.max_heavies
        if max_heavies is not None and num_heavies > max_heavies:
            return "too many heavy atoms"

        max_rotatable_bonds = self.max_rotatable_bonds
        if max_rotatable_bonds is not None:
            # The maxMatches (in RDKit at least up to early 2016) specifies
            # the number of non-unique matches. After those N are found, they
            # are filtered to remove duplicates. It's possible that up to
            # 1/2 of the bonds will be removed, so the limit must be twice
            # the max requested, plus 1
            matches = mol.GetSubstructMatches(self.rotatable_pattern, uniquify=True,
                                              maxMatches=max_rotatable_bonds*2+1)
            if len(matches) > max_rotatable_bonds:
                return "too many rotatable bonds"

        return None

    def get_cut_atom_pairs(self, mol):
        return mol.GetSubstructMatches(self.cut_pattern, uniquify=True)
                    
    def get_cut_lists(self, mol):
        atom_pairs = self.get_cut_atom_pairs(mol)

        # Generate all of the atoms lists with 1, 2, or 3 cuts
        cut_lists = []
        for i, first_pair in enumerate(atom_pairs):
            cut_lists.append( [first_pair] )
            
            if self.num_cuts >= 2:
                for j, second_pair in enumerate(atom_pairs[i+1:], i+1):
                    cut_lists.append( [first_pair, second_pair] )
                    
                    if self.num_cuts >= 3:
                        for k, third_pair in enumerate(atom_pairs[j+1:], j+1):
                            cut_lists.append( [first_pair, second_pair, third_pair] )
        return cut_lists

class FragmentValueError(ValueError):
    def __init__(self, name, value, reason):
        self.name = name
        self.value = value
        self.reason = reason

    def __str__(self):
        return "Error with %s (%r): %s" % (self.name, self.value, self.reason)

    def __repr__(self):
        return "FragmentValueError(%r, %r, %r)" % (
            self.name, self.value, self.reason)

    
def get_fragment_filter(fragment_options):
    options = fragment_options
    
    def call(parse, name):
        value = getattr(options, name)
        try:
            return parse(value)
        except ValueError as err:
            raise FragmentValueError(name, value, str(err))
    
    max_heavies = options.max_heavies
    max_rotatable_bonds = options.max_rotatable_bonds
    rotatable_pattern = call(parse_rotatable_smarts, "rotatable_smarts")
    cut_pattern = call(parse_cut_smarts, "cut_smarts")
    num_cuts = call(parse_num_cuts, "num_cuts")
    method = call(parse_method, "method")
    try:
        salt_remover = call(parse_salt_remover, "salt_remover")
    except IOError as err:
        raise FragmentValueError("salt_remover", options.salt_remover,
                                 "Cannot open salt file: %s" % (err,))
    
    return FragmentFilter(
        max_heavies = max_heavies,
        max_rotatable_bonds = max_rotatable_bonds,
        rotatable_pattern = rotatable_pattern,
        salt_remover = salt_remover,
        cut_pattern = cut_pattern,
        num_cuts = num_cuts,
        method = method,
        options = fragment_options,
        )
    

###


class ParsedSmilesRecord(object):
    __slots__ = ("id", "smiles", "mol", "normalized_mol", "normalized_smiles", "num_normalized_heavies")
    def __init__(self, id, smiles, mol, normalized_mol, normalized_smiles, num_normalized_heavies):
        self.id = id
        self.smiles = smiles
        self.mol = mol
        self.normalized_mol = normalized_mol
        self.normalized_smiles = normalized_smiles
        self.num_normalized_heavies = num_normalized_heavies


def parse_record(id, smiles, fragment_filter):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "invalid smiles", ParsedSmilesRecord(id, smiles, mol, None, None, 0)

    errmsg, normalized_mol = fragment_filter.normalize(mol)
    normalized_smiles = Chem.MolToSmiles(normalized_mol, isomericSmiles=True)
    if errmsg is None:
        if "." in normalized_smiles:
            errmsg = "multiple fragments"
    num_normalized_heavies = fragment_algorithm.count_num_heavies(mol)

    record = ParsedSmilesRecord(id, smiles, mol, normalized_mol, normalized_smiles, num_normalized_heavies)
    if errmsg is not None:
        return errmsg, record

    errmsg = fragment_filter.apply_filters(normalized_mol)
    if errmsg is not None:
        return errmsg, record

    return None, record


###


# Adapater to emulate the multiprocessing pool without using any threads/processes
class SingleProcessPool(object):
    def apply_async(self, f, args):
        return SyncResult(f, args)
    def terminate(self):
        pass
    def join(self):
        pass
    def close(self):
        pass

class SyncResult(object):
    def __init__(self, f, args):
        self.f = f
        self.args = args
        
    def get(self):
        # Don't compute until requested
        # Note: if you .get() twice, I'll compute twice. That's all that
        # this code needs.
        return self.f(*self.args)

# Make it so that ^C works
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def create_pool(num_jobs):
    if num_jobs > 1:
        pool = multiprocessing.Pool(num_jobs, init_worker)
    else:
        pool = SingleProcessPool()
    return pool
    
####


def _as_list(method, normalized_mol, fragment_filter, num_normalized_heavies):
    return list(method(normalized_mol, fragment_filter, num_normalized_heavies))

def make_hydrogen_fragment_record(id, input_smiles, fragment_filter):
    errmsg, record = parse_record(id, input_smiles, fragment_filter)
    if errmsg:
        return fragment_types.FragmentErrorRecord(id, input_smiles, errmsg)

    fragments = fragment_algorithm.fragment_molecule_on_explicit_hydrogens(input_smiles)
    return fragment_types.FragmentRecord(
        id, input_smiles,
        record.num_normalized_heavies, record.normalized_smiles,
        fragments)
    


def make_fragment_records(smiles_reader, fragment_filter, cache=None, pool=None, reporter=None):
    reporter = reporters.get_reporter(reporter)
    jobs = []

    if pool is None:
        pool = SingleProcessPool()


    # There are two phases:
    #   1) establish what needs to be fragmented vs. what is available from
    #       cache or could not be parsed
    #   2) fragment the unfragmented
    for recno, terms in reporter.progress(
            enumerate(smiles_reader), "Preparing record"):
        input_smiles = terms[0]
        id = terms[1]
        where = smiles_reader.location.where()

        # If the fragment information is available from cache then use it
        if cache is not None:
            record = cache.get(id)
            if record is not None:
                if record.input_smiles == input_smiles:
                    jobs.append( (id, input_smiles, where, None, record) )
                    continue

        # If I can't parse it then record the error messages
        errmsg, record = parse_record(id, input_smiles, fragment_filter)
        if errmsg:
            result = fragment_types.FragmentErrorRecord(id, input_smiles, errmsg)
            jobs.append( (id, input_smiles, where, None, result) )
            continue

        # Submit it as something to work on
        args = (fragment_filter.method, record.normalized_mol, fragment_filter,
                record.num_normalized_heavies)
        result = pool.apply_async(_as_list, args)

        jobs.append( (id, input_smiles, where, record, result) )

    # I'll a bit cautious. I'll process the jobs in order, yield
    # the result, then throw it away. This keeps the job list from
    # being filled with completed results.
    def pop_iter(jobs):
        while jobs:
            yield jobs.pop(0)
    with reporter.progress(pop_iter(jobs),
                           "Fragmented record", len(jobs)) as job_iter:
        for (id, input_smiles, where, record, result) in job_iter:
            if record is None:
                # use a pre-computed result (from cache, or an error record)
                yield result
                continue

            try:
                fragments = result.get()
            except RuntimeError as err:
                # Some sort of RDKit failure.
                reporter.update("")
                print("Skipping %s: %s" % (err, where), file=sys.stderr)
                continue

            except Exception:
                # Something unexpected happened.
                # Give some idea of what failed.
                reporter.update("")
                print("Failure:", where, file=sys.stderr)
                raise

            yield fragment_types.FragmentRecord(
                id, input_smiles,
                record.num_normalized_heavies, record.normalized_smiles,
                fragments)


class SingleSmilesReader(object):
    def __init__(self, smiles, id="query"):
        self.id = id
        self.smiles = smiles
    def __iter__(self):
        yield (self.smiles, self.id)
    location = fileio.Location.from_source("<string>")

def make_fragment_record_from_smiles(smiles, fragment_filter):
    reader = SingleSmilesReader(smiles)
    records = make_fragment_records(reader, fragment_filter)
    for record in records:
        return record
    raise AssertionError("how can there not be any records?")
        

###
def _name_to_command_line(s):
    return "--" + s.replace("_", "-")

def cannot_combine_cache_and_fragment_args(parser, specified_args):
    names = sorted(_name_to_command_line(name) for name in specified_args)
    if len(names) == 1:
        parser.error("Cannot combine %r with --cache" % (names[0],))
    else:
        left = ", ".join(map(repr, names[:-1]))
        parser.error("Cannot combine %s or %r with --cache" % (left, names[-1]))


def fragment_command(parser, args):
    # How do I report errors?
    reporter = command_support.get_reporter(args.quiet)
    
    # Get the command-line fragment arguments
    specified_args, options = get_fragment_options_from_args(parser, args)
    
    # Use a cache?
    cache = None
    if args.cache is not None:
        if specified_args:
            # Must use the arguments from the cache.
            # Cannot specify command-line arguments.
            cannot_combine_cache_and_fragment_args(parser, specified_args)
            raise AssertionError("Should not get here")
        
        try:
            cache = fragment_io.load_cache(args.cache, reporter)
        except IOError as err:
             parser.error("Cannot open cache: %s" % (err,))
        except ValueError as err:
             parser.error("Problem loading cache: %s" % (err,))

        try:
            fragment_filter = get_fragment_filter(cache.options)
        except FragmentValueError as err:
            parser.error("Error in cache option %r (%r) from %r: %s" % (
                err.name, err.value, args.cache, err.reason))
    else:
        try:
            fragment_filter = get_fragment_filter(options)
        except FragmentValueError as err:
            parser.error("Error in command-line option %r (%r): %s" % (
                _name_to_command_line(err.name), err.value, err.reason))

    pool = create_pool(args.num_jobs)
    # We need to re-enable catching ^C in order to kill off the pool gracefully.
    # (^C was disabled in 'mmpdb')
    signal.signal(signal.SIGINT, signal.default_int_handler)
    
    try:
        try:
            with fileio.read_smiles_file(args.structure_filename, args.format,
                                         args.delimiter, args.has_header) as reader:
                with fragment_io.open_fragment_writer(filename = args.output,
                                                      options = fragment_filter.options,
                                                      format_hint = args.out) as writer:
                    records = make_fragment_records(reader, fragment_filter, cache,
                                                    pool=pool, reporter=reporter)
                    writer.write_records(records)

        except FileFormatError as err:
            sys.stderr.write("Cannot parse input file: %s\n" % (err,))
            raise SystemExit(1)

    except KeyboardInterrupt:
        reporter.update("Shutting down process pool")
        pool.terminate()
        pool.join()
        reporter.update("")
        raise SystemExit(-1)
    else:
        reporter.update("Closing process pool")
        pool.close()
        pool.join()
        reporter.update("")
    
        ## with fileio.open_output(args.output, args.out) as outfile:
        ##     fragment_io.write_fragment_records(outfile, records, )

def smifrag_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)
    specified, fragment_options = get_fragment_options_from_args(parser, args)
    fragment_filter = get_fragment_filter(fragment_options)
    record = make_fragment_record_from_smiles(args.smiles, fragment_filter)
    if record.errmsg:
        parser.error("Cannot parse --smiles: %s" % (record.errmsg,))
        
    #writer = fragment_io.FragInfoWriter(None, sys.stdout, None)
    #writer.write_records([record])

    columns = [["#cuts"], ["enum.label"],
               ["#heavies"], ["symm.class"], ["smiles"],
               ["order"],
               ["#heavies"], ["symm.class"], ["smiles"],
               ["with-H"]]
    styles = ["center", "center",
              "right", "center", "left", "center",
              "right", "center", "left", "left"]

    has_rows = False
    for frag in record.fragments:
        has_rows = True
        items = [str(frag.num_cuts),
                 frag.enumeration_label,
                 str(frag.variable_num_heavies),
                 frag.variable_symmetry_class,
                 frag.variable_smiles,
                 frag.attachment_order,
                 str(frag.constant_num_heavies),
                 frag.constant_symmetry_class,
                 frag.constant_smiles,
                 frag.constant_with_H_smiles or "-",
                 ]
                     
        for (item, column) in zip(items, columns):
            column.append(str(item))

    sizes = []
    for style, column in zip(styles, columns):
        column_width = max(map(len, column))
        column[0] = column[0].ljust(column_width)
        sizes.append(column_width)
        if len(column) == 1:
            continue
        
        data_width = max(map(len, column[1:]))
        
        if style == "center":
            column[1:] = [s.rjust(data_width).center(column_width) for s in column[1:]]
        elif style == "left":
            column[1:] = [s.ljust(data_width) for s in column[1:]]
        ## elif style == "right-10":
        ##     spacer = " "*10
        ##     column[1:] = [spacer + s.rjust(data_width-10).center(column_width-10) for s in column[1:]]
        elif style == "right":
            column[1:] = [s.rjust(data_width).center(column_width) for s in column[1:]]
        else:
            raise AssertionError(style)

    first_line = (
        " " * sizes[0] + "   " +
        " " * sizes[1] + " |-" +
        "  variable  ".center(sizes[2] + sizes[3] + sizes[4] + 6, "-") + "-| " +
        " " * sizes[5] + " |-" +
        "  constant  ".center(sizes[6] + sizes[7] + sizes[8] + sizes[9] + 9, "-")
        )

    print(first_line)
    for lineno, fields in enumerate(zip(*columns)):
        print(*fields, sep=" | ")
        if lineno == 0:
            print(*["-"*len(s) for s in fields], sep = "-+-")
        
