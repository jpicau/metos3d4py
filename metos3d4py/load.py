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

        # sequential run
        text = ""
        if size == 1:
            text = text + "sequential run, 1 process"
            return text

        # parallel run
        if rank == 0:
            header = ""
            header = header + "parallel run, {} processes\n".format(size)
            header = header + "process   profiles                    vector"

            opt = self.opt
            maxdiff = self.maxdiff
            
            footer = ""
            footer = footer + "optimal vector length: {:.1f}\n".format(opt)
            footer = footer + "   max relative error: {:2.6f}\n".format(maxdiff/opt)
            footer = footer + "   max absolute error: {:.1f}".format(maxdiff)

        text = text + "{:3d}/{:<3d}   ".format(rank, size)
        text = text + "start, count: {:5d} {:5d}   ".format(self.npprev, self.nploc)
        text = text + "start, count: {:7d} {:7d}   ".format(self.nvprev, self.nvloc)

        if size > 1:
            msgs = []
            if rank == 0:
                msgs.append(header)
                msgs.append(text)
                for i in range(size-1):
                    # receive
                    req = comm.irecv(source=i+1, tag=i+1)
                    msgs.append(req.wait())
                msgs.append(footer)
                return "\n".join(msgs)
            else:
                # send
                req = comm.isend(text, dest=0, tag=rank)
                req.wait()
    
        return text

# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):
        
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
    
        # profiles to processes
        starts = npi.cumsum() - npi
        weights = starts + 0.5*npi
        ranks = numpy.floor((weights/nv)*size)

        # profiles per process
        punique, pindex, pcount = numpy.unique(ranks, return_index = True, return_counts = True)
        pprev = pcount.cumsum() - pcount

        # vector shares
        vcount = numpy.array([npi[i:i+c].sum() for i, c in zip(pindex, pcount)])
        vprev = vcount.cumsum() - vcount
        
        # my rank
        self.nploc = pcount[rank]
        self.npprev = pprev[rank]
        self.nvloc = vcount[rank]
        self.nvprev = vprev[rank]
        
        
        # index pointers for bgc
#        numpy.append(0, numpy.cumsum(nnzrow)).astype("i4")
        self.npiloc = npi[self.npprev:self.npprev+self.nploc]
        self.indptr = numpy.append(0, numpy.cumsum(npi))
        self.indptrloc = numpy.append(0, numpy.cumsum(self.npiloc))
#        self.np_indptr = numpy.append(0, numpy.cumsum(npi[self.npprev:self.npprev+self.nploc]))
#        self.

        # compute max diff to optimal
        self.opt = opt = float(nv)/size
        self.maxdiff = maxdiff = numpy.amax(numpy.abs(numpy.array(vcount) - opt))

        # debug
        util.debug(m3d, self, self, level=1)








