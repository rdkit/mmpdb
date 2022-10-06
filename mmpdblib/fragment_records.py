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

from . import fragment_types
from . import fileio


def _as_list(method, normalized_mol, fragment_filter, num_normalized_heavies):
    return list(method(normalized_mol, fragment_filter, num_normalized_heavies))


###


class ParsedSmilesRecord(object):
    __slots__ = (
        "id",
        "smiles",
        "mol",
        "normalized_mol",
        "normalized_smiles",
        "num_normalized_heavies",
    )

    def __init__(self, id, smiles, mol, normalized_mol, normalized_smiles, num_normalized_heavies):
        self.id = id
        self.smiles = smiles
        self.mol = mol
        self.normalized_mol = normalized_mol
        self.normalized_smiles = normalized_smiles
        self.num_normalized_heavies = num_normalized_heavies


def parse_record(id, smiles, fragment_filter):
    from rdkit import Chem
    from . import fragment_algorithm

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "invalid smiles", ParsedSmilesRecord(id, smiles, mol, None, None, 0)

    errmsg, normalized_mol = fragment_filter.normalize(mol)
    normalized_smiles = Chem.MolToSmiles(normalized_mol, isomericSmiles=True)
    if errmsg is None:
        if "." in normalized_smiles:
            errmsg = "multiple fragments"
    num_normalized_heavies = mol.GetNumHeavyAtoms()

    record = ParsedSmilesRecord(id, smiles, mol, normalized_mol, normalized_smiles, num_normalized_heavies)
    if errmsg is not None:
        return errmsg, record

    errmsg = fragment_filter.apply_filters(normalized_mol)
    if errmsg is not None:
        return errmsg, record

    return None, record


def make_hydrogen_fragment_record(id, input_smiles, fragment_filter):
    from . import fragment_algorithm

    errmsg, record = parse_record(id, input_smiles, fragment_filter)
    if errmsg:
        return fragment_types.FragmentErrorRecord(id, input_smiles, errmsg)

    fragments = fragment_algorithm.fragment_molecule_on_explicit_hydrogens(input_smiles)
    return fragment_types.FragmentRecord(
        id,
        input_smiles,
        record.num_normalized_heavies,
        record.normalized_smiles,
        fragments,
    )


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


def make_fragment_records(smiles_reader, fragment_filter, cache=None, pool=None, reporter=None):
    jobs = []

    if pool is None:
        pool = SingleProcessPool()

    # There are two phases:
    #   1) establish what needs to be fragmented vs. what is available from
    #       cache or could not be parsed
    #   2) fragment the unfragmented
    for recno, terms in reporter.progress(enumerate(smiles_reader), "Preparing record"):
        input_smiles = terms[0]
        id = terms[1]
        where = smiles_reader.location.where()

        # If the fragment information is available from cache then use it
        if cache is not None:
            record = cache.get(id)
            if record is not None:
                if record.input_smiles == input_smiles:
                    jobs.append((id, input_smiles, where, None, record))
                    continue

        # If I can't parse it then record the error messages
        errmsg, record = parse_record(id, input_smiles, fragment_filter)
        if errmsg:
            result = fragment_types.FragmentErrorRecord(id, input_smiles, errmsg)
            jobs.append((id, input_smiles, where, None, result))
            continue

        # Submit it as something to work on
        args = (
            fragment_filter.method,
            record.normalized_mol,
            fragment_filter,
            record.num_normalized_heavies,
        )
        result = pool.apply_async(_as_list, args)  # fragment_filter.method calls the actual fragmentation algorithm

        jobs.append((id, input_smiles, where, record, result))

    # I'll a bit cautious. I'll process the jobs in order, yield
    # the result, then throw it away. This keeps the job list from
    # being filled with completed results.
    def pop_iter(jobs):
        while jobs:
            yield jobs.pop(0)

    with reporter.progress(pop_iter(jobs), "Fragmented record", len(jobs)) as job_iter:
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

            except fragment_types.FragmentationFailure as err:
                yield fragment_types.FragmentErrorRecord(
                    id,
                    input_smiles,
                    str(err),
                    )
                continue

            except Exception:
                # Something unexpected happened.
                # Give some idea of what failed.
                reporter.update("")
                print("Failure:", where, file=sys.stderr)
                raise

            yield fragment_types.FragmentRecord(
                id,
                input_smiles,
                record.num_normalized_heavies,
                record.normalized_smiles,
                fragments,
            )


########


class SingleSmilesReader(object):
    def __init__(self, smiles, id="query"):
        self.id = id
        self.smiles = smiles

    def __iter__(self):
        yield (self.smiles, self.id)

    location = fileio.Location.from_source("<string>")


def make_fragment_record_from_smiles(smiles, fragment_filter, reporter=None):
    from . import reporters

    reporter = reporters.get_reporter(reporter)

    reader = SingleSmilesReader(smiles)
    records = make_fragment_records(reader, fragment_filter, reporter=reporter)
    for record in records:
        return record
    raise AssertionError("how can there not be any records?")
