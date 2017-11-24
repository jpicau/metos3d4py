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
        
#        if size > 1:
#            util.debug(m3d, self, "parallel run, {} processes".format(size), level=1)

        text = text + "rank: {:3d}/{}   ".format(rank, size)
        text = text + "npprev, nploc: {:5d} {:5d}   ".format(self.npprev, self.nploc)
        text = text + "nvprev, nvloc: {:7d} {:7d}   ".format(self.nvprev, self.nvloc)
        
        if rank == 0:
            text = "distribution: \n" + text
        
        if size > 1:
            msgs = []
            if rank == 0:
#                text = "Load:\n" + text
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
        
        # compute max diff to optimal
        self.opt = opt = float(nv)/size
        self.maxdiff = maxdiff = numpy.amax(numpy.abs(numpy.array(vcount) - opt))

#        print("{} {}, {} {}".format(self.nvloc, self.nvprev, , ))

#        # store self in m3d
#        m3d.load = self

#        # debug
#        if size > 1:
#            util.debug(self, self, "parallel run, {} processes".format(size), level=1)
#        else:
#            util.debug(self, self, "sequential run, {} process".format(size), level=1)
#        util.debug(self, self, "config file: {}".format(argv[1]), level=1)

#        if size > 1:
#            util.debug(m3d, self, "parallel run, {} processes".format(size), level=1)
#
#            # compute max diff to optimal
#            self.opt = opt = float(nv)/size
#            self.maxdiff = maxdiff = numpy.amax(numpy.abs(numpy.array(vcount) - opt))
#            util.debug(
#                       m3d,
#                       self,
#                       "optimum: {:.1f}, max error: relative: {:2.6f}, absolute: {}"
#                       .format(opt, maxdiff/opt, int(maxdiff)), level=1)
#        else:
#            util.debug(m3d, self, "sequential run, {} process".format(size), level=1)

        # debug
        util.debug(m3d, self, self, level=1)



# ----------------------------------------------------------------------------------------




