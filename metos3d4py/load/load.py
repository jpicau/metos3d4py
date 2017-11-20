#
# Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
# Copyright (C) 2017  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy

"""
    Load submodule
    ==============
    
    Attributes:
        nploc
        npprev
        nvloc
        nvprev
        
    """
class Load:
    def init(self, comm, grid):
        size = comm.size
        rank = comm.rank
        
        nv = grid.nv
        np = grid.np
        npi = grid.npi
    
        # assign profiles to processes
        starts = npi.cumsum() - npi
        weights = starts + 0.5*npi
        ranks = numpy.floor((weights/nv)*size)

        rank_unique, rank_index, rank_count = numpy.unique(ranks, return_index = True, return_counts = True)
        rank_prev = rank_count.cumsum() - rank_count

        nploc = rank_count[rank]
        npprev = rank_prev[rank]

        nvloc = npi[npprev:npprev+nploc].sum()
        nvprev = npi[0:npprev].sum()

        self.nploc = nploc
        self.npprev = npprev
        self.nvloc = nvloc
        self.nvprev = nvprev

    def __str__(self):
        return "Load:\n  npprev, nploc: {} {}\n  nvprev, nvloc: {} {}".format(self.npprev, self.nploc, self.nvprev, self.nvloc)


