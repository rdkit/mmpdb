"index configuration"
# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
# Copyright (c) 2021, F. Hoffmann-La Roche Ltd.
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

# These are shared by a few different modules and placed into its own
# module to reduce inter-dependencies between them.


class IndexOptions(object):
    __slots__ = (
        "min_variable_heavies",
        "max_variable_heavies",
        "min_variable_ratio",
        "max_variable_ratio",
        "min_radius",
        "max_radius",
        "symmetric",
        "max_heavies_transf",
        "max_frac_trans",
        "smallest_transformation_only",
    )

    def __init__(
        self,
        min_variable_heavies = None,
        max_variable_heavies = None,
        min_variable_ratio = None,
        max_variable_ratio = None,
        max_heavies_transf = None,
        max_frac_trans = None,
        min_radius = 0,
        max_radius = 5,
        symmetric = False,
        smallest_transformation_only = False,
    ):

        assert min_variable_heavies is None or min_variable_heavies >= 0, min_variable_heavies
        self.min_variable_heavies = min_variable_heavies

        assert (
            (max_variable_heavies is None)
            or (min_variable_heavies is None and max_variable_heavies >= 0)
            or (min_variable_heavies <= max_variable_heavies)
        ), max_variable_heavies
        self.max_variable_heavies = max_variable_heavies

        assert min_variable_ratio is None or 0.0 <= min_variable_ratio <= 1.0, min_variable_ratio
        self.min_variable_ratio = min_variable_ratio

        assert (
            (max_variable_ratio is None)
            or (min_variable_ratio is None and max_variable_ratio <= 1.0)
            or (min_variable_ratio <= max_variable_ratio <= 1.0)
        )
        self.max_variable_ratio = max_variable_ratio

        assert max_heavies_transf is None or max_heavies_transf >= 0, max_heavies_transf
        self.max_heavies_transf = max_heavies_transf

        assert max_frac_trans is None or max_frac_trans >= 0, max_heavies_transf
        self.max_frac_trans = max_frac_trans

        assert min_radius <= max_radius, (min_radius, max_radius)
        
        assert min_radius >= 0, min_radius
        self.min_radius = min_radius

        assert max_radius >= 0, max_radius
        self.max_radius = max_radius

        assert isinstance(symmetric, bool)
        self.symmetric = symmetric

        assert isinstance(smallest_transformation_only, bool)
        self.smallest_transformation_only = smallest_transformation_only

    def to_dict(self):
        d = {}
        for name in IndexOptions.__slots__:
            value = getattr(self, name)
            if value is not None:
                d[name] = value
        return d

    def get_fragment_filter(self):
        from . import index_algorithm

        filters = []
        if self.min_variable_heavies is not None:
            filters.append(index_algorithm.MinVariableHeaviesFilter(self.min_variable_heavies))
        if self.max_variable_heavies is not None:
            filters.append(index_algorithm.MaxVariableHeaviesFilter(self.max_variable_heavies))
        if self.min_variable_ratio is not None:
            filters.append(index_algorithm.MinVariableRatioFilter(self.min_variable_ratio))
        if self.max_variable_ratio is not None:
            filters.append(index_algorithm.MaxVariableRatioFilter(self.max_variable_ratio))

        if not filters:
            # It's easier to have 0 filters than to make a special do-nothing filter.
            return index_algorithm.MultipleFilters([])
        elif len(filters) == 1:
            return filters[0]
        else:
            return index_algorithm.MultipleFilters(filters)
