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
        optimal load balancing,
        for given grid, mask, profiles, vector,
        between processes,
        
        Attributes:
        nploc
        npprev
        nvloc
        nvprev
        
    """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        """
            Collect messages from each process in an ordered/synchronized manner.
            Use the internal mpi4py communicator.
            Rank 0 receives, all other process send.
        """
        comm = self.comm.tompi4py()
        size = comm.size
        rank = comm.rank

        text = ""
        text = text + "  rank: {:3d}/{}   ".format(rank, size)
        text = text + "  npprev, nploc: {:5d} {:5d}   ".format(self.npprev, self.nploc)
        text = text + "  nvprev, nvloc: {:7d} {:7d}   ".format(self.nvprev, self.nvloc)
        if size > 1:
            msgs = []
            if rank == 0:
                text = "Load:\n" + text
                msgs.append(text)
                for i in range(size-1):
                    # receive
                    req = comm.irecv(source=i+1, tag=i+1)
                    msgs.append(req.wait())
                return "\n".join(msgs)
            else:
                # send
                req = comm.isend(text, dest=0, tag=rank)
                req.wait()
    
        return text

# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        
        comm = m3d.comm
        grid = m3d.grid

        # store own copy of comm for __str__
        self.comm = m3d.comm

        self.size = size = comm.size
        self.rank = rank = comm.rank

        # vector length, profile count, profile depths
        nv = grid.nv
        np = grid.np
        npi = grid.npi
    
        # assign profiles to processes
        starts = npi.cumsum() - npi
        weights = starts + 0.5*npi
        ranks = numpy.floor((weights/nv)*size)

        # count profiles per process
        punique, pindex, pcount = numpy.unique(ranks, return_index = True, return_counts = True)
        pprev = pcount.cumsum() - pcount

        # store own in self
        self.nploc = nploc = pcount[rank]
        self.npprev = npprev = pprev[rank]
        self.nvloc = npi[npprev:npprev+nploc].sum()
        self.nvprev = npi[0:npprev].sum()

        # store self in m3d
        m3d.load = self

        if size > 1:
            # count vector shares
            vcount = [npi[i:i+c].sum() for i, c in zip(pindex, pcount)]
            # compute max diff to optimal
            opt = float(nv)/size
            maxdiff = numpy.amax(numpy.abs(numpy.array(vcount) - opt))
            util.debug(
                       m3d,
                       self,
                       "Optimum: {:.1f}, max error: relative: {:2.6f}, absolute: {}"
                       .format(opt, maxdiff/opt, int(maxdiff)), level=1)

# ----------------------------------------------------------------------------------------




