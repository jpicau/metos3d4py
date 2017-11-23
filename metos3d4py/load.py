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
from metos3d4py import util

# ----------------------------------------------------------------------------------------
class Load:
    """
        Load submodule
        ==============
        
        Attributes:
        nploc
        npprev
        nvloc
        nvprev
        
    """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        size = self.size
        rank = self.rank
        
        text = ""
        if rank == 0:
            text = text + "Load:\n"
        text = text + "  rank: {:3d}/{}   ".format(rank, size)
        text = text + "  npprev, nploc: {:5d} {:5d}   ".format(self.npprev, self.nploc)
        text = text + "  nvprev, nvloc: {:7d} {:7d}   ".format(self.nvprev, self.nvloc)
        return text

# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        
        comm = m3d.comm
        grid = m3d.grid
        util.debug(m3d, "Load init: {}".format("..."), level=1)

        # store own copy of rank and size for __str__
        self.size = size = comm.size
        self.rank = rank = comm.rank
        
        # ...
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

        # store in self
        self.nploc = nploc
        self.npprev = npprev
        self.nvloc = nvloc
        self.nvprev = nvprev
    
        # store self in m3d
        m3d.load = self



