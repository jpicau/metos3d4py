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

"""
    Metos3D for Python
    ==================
    
    This is a Python version of Metos3D.
    
    used packages:
        h5py        hdf5 for python, netdcf 4 files, binary, parallel i/o
        petcs4py    parallel computation, data types, operations, solvers, ...
        yaml        option format, compact, human readable and writable
"""

from petsc4py import PETSc

from metos3d4py.version     import VERSION
from metos3d4py.grid        import Grid
from metos3d4py.load        import Load
from metos3d4py.tracer      import Tracer
from metos3d4py.bgc         import BGC
from metos3d4py.tmm         import TMM
from metos3d4py.time        import Time
from metos3d4py.solver      import Solver
from metos3d4py             import util

# ----------------------------------------------------------------------------------------
class Metos3D:
    """
        Metos3D class
        
        serves as the Metos3D data type, i.e.
        stores the metos3d simulation context,

        Main attributes:
            grid:
            load:
            tracer:
            bgc:
            tmm:
            time:
            solve:

        Methods:
            init():     initialize metos3d context from a yaml configuration file
            run():      perform a simulation
            final():    finalize the context
    
    """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        return "metos3d version {} ".format(VERSION)

# ----------------------------------------------------------------------------------------
    def __init__(self, argv):
        """
            check computation model ...
            print version,
            print number of processes,
            attributes ...
            creates instances ...
            initialize metos3d context,

            Parameters:
                argv:   list with command arguments, file path of yaml conf file

        """
        
        # version and runtime context
        self.version    = version   = VERSION
        self.comm       = comm      = PETSc.COMM_WORLD
        self.size       = size      = comm.size
        self.rank       = rank      = comm.rank
        
        # config
        self.config = util.get_config_from_yaml_file(self, argv)

        # debug level
        self.debug = util.get_key(self, self, self.config, "Debug", int)
        util.debug(self, self, self, level=1)

        # init
        self.grid = Grid(self)
        self.load = Load(self)
#        self.load = load = Load(self)
#        util.debug(self, load, load, level=1)

#, load, tracer, bgc, tmm, time, solver = \
#            Grid(), Load(), Tracer(), BGC(), TMM(), Time(), Solver()

#        # init
#        grid.init(self)
#        load.init(self)
#        tracer.init(self)
#        bgc.init(self)
#        tmm.init(self)
#        time.init(self)
#        solver.init(self)

        # debug



# ----------------------------------------------------------------------------------------
    def run(self):
        """
            runs a simulation,
            loop over model years, spin up
            perform time steps,
            interpolate, boundary and domain data, matrices,
            evaluate bgc model,
            apply transport,
            
            u = [u1, ..., um]
            bj = [[]]
            dj = [[]]
            yj = [yj1, ..., yjny]
            qj = [qj1, ..., qjny]
            Aimpj = diag([Aimpj, ..., Aimpj])
            Aexpj = diag([Aexpj, ..., Aexpj])
            
            yjp1 = Aimpj * (Aexpj * yj + qj(dt, tj, yj, u, bj, dj))
            
            yjp1 = [yjp11, ..., yjp1ny]
            
        """
        
        pass
        # compute solution
#        self.solver.solve(self)

        # store as hdf5 file
#        self.tracer.write(self)

#        comm = self.comm
#        t0, dt, nt = self.time.get()
#        ny, y0, yj, yexpj, qj = self.bgc.tracer.get()
##        ny = self.bgc.ny
##        y0 = self.bgc.y0
##        yj = self.bgc.yj
##        yexpj = self.bgc.yexpj
##        qj = self.bgc.qj
#
#        nu = self.bgc.nu
#        u = self.bgc.u
#
#        bgc = self.bgc.bgc
#
#        nb = self.bgc.nb
#        nbi = self.bgc.nbi
#        b = self.bgc.b
#        bj = self.bgc.bj
#
#        if hasattr(self.bgc, "nd"):
#            nd = self.bgc.nd
#            ndi = self.bgc.ndi
#            d = self.bgc.d
#            dj = self.bgc.dj
#        else:
#            nd = 0
#            ndi = []
#            d = []
#            dj = []

#        nl = self.solver.nl
#        yl = self.solver.yl
#
#        nexp = self.tmm.nexp
#        nimp = self.tmm.nimp
##        Aexp, Aimp
#
#        # init
#        for i in range(ny):
#            yl[i] = y0[i]
#
#        # spin up
#        for l in range(nl):
#
#            # store state at begin of model year
#            for i in range(ny):
#                yj[i] = yl[i]
#
#            # time step
#            for j in range(nt):
#                tj = t0 + j*dt
#                util._print(comm, "Time step: {}/{} {}".format(j, nt, tj))
#
#                # interpolate
#                for ib in range(nb):
#                    alpha, ialpha, beta, ibeta = _interpolate(nbi[ib], tj)
#                    bj[ib] = alpha * b[ib][ialpha] + beta * b[ib][ibeta]
#
#                for id in range(nd):
#                    alpha, ialpha, beta, ibeta = _interpolate(ndi[id], tj)
#                    dj[id] = alpha * d[id][ialpha] + beta * d[id][ibeta]
#
#                alpha, ialpha, beta, ibeta = _interpolate(nexp, tj)
#                Aexpj = alpha * Aexp[ialpha] + beta * Aexp[ibeta]
#
#                alpha, ialpha, beta, ibeta = _interpolate(nimp, tj)
#                Aimpj = alpha * Aimpj[ialpha] + beta * Aimpj[ibeta]
#
#                # bgc
#                bgc(self, dt, qj, tj, yj, u, bj, dj)
#
#                # tmm
#                for i in range(ny):
#                    yexpj[i] = Aexpj * yj[i] + qj[i]
#                    yj[i]    = Aimpj * yexpj[i]
#
#                # error
##                norm2 = norm2(yl, yj)

# ----------------------------------------------------------------------------------------
    def final(self):
        """
            final()

            finalize metos3d context,
        """
        
        # if output
        # store to disc
        
        pass
















