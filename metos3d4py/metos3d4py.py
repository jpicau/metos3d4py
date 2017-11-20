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

from metos3d4py             import VERSION
from metos3d4py.conf        import Conf
from metos3d4py.grid        import Grid
from metos3d4py.load        import Load
from metos3d4py.bgc         import Bgc
from metos3d4py.transport   import Transport
from metos3d4py.util.util   import _print_message, _print_message_synch

class Metos3D:
    """
        Metos3D class
        
        serves as the Metos3D data type, i.e.
        stores the metos3d simulation context,

        Attributes:
            version
            
            comm
            size
            rank
            
            conf:
            grid:
            load:
            bgc:

        Methods:
            init():     initialize metos3d context from a yaml configuration file
            run():      perform a simulation
            final():    finalize the context
    
    """

# ----------------------------------------------------------------------------------------

    def final(self):
        """
            final()
            
            finalize metos3d context,
            """
        pass

# ----------------------------------------------------------------------------------------

    def init(self, argv):
        """
            init()
            
            Parameters:
                argv:   list with command arguments
            
            print version,
            print number of processes,

            initialize metos3d context,
                conf,
                grid,
                load,
                bgc,
                
            """
        
        self.version    = version   = VERSION
        self.comm       = comm      = PETSc.COMM_WORLD
        self.size       = size      = comm.size
        self.rank       = rank      = comm.rank
                    
        _print_message(comm, "metos3d version {} ".format(version))
        if size > 1:
            _print_message(comm, "parallel run, {} processes".format(size))
        else:
            _print_message(comm, "sequential run, {} process".format(size))

        self.conf       = conf      = Conf()
        self.grid       = grid      = Grid()
        self.load       = load      = Load()
        self.bgc        = bgc       = Bgc()
        self.transport  = transport = Transport()
        
        conf.init(comm, argv)
        grid.init(comm, conf)
        load.init(comm, grid)
        bgc.init(comm, conf, grid, load)
        transport.init(comm, conf, grid, load)
        
        # debug
#        _print_message(comm, conf)
#        _print_message(comm, grid)
#        _print_message_synch(comm, load)
#        _print_message(comm, bgc)


## transport ------------------------------------------------------------------------------------
###self.Aexp, nexp
###    self.Aexpj
##        Aexp, nexp = _get_transport_explicit_conf()
##        Aexpj = Aexp.duplicate()
###self.Aimp, nimp
###    self.Aimpj
##        Aimp, nimp = _get_transport_implicit_conf()
##        Aimpj = Aimp.duplicate()
#        transport_list = self.conf.dict["Transport matrix"]["Name, Count, Type, File"]
#        transport_patj = self.conf.dict["Transport matrix"]["Path"]
#
#
## time step
##        t0, nt, dt = _get_time_step_conf()
#        self.t0 = self.conf.dict["Time step"]["Start"]
#        self.nt = self.conf.dict["Time step"]["Count"]
#        self.dt = self.conf.dict["Time step"]["Step"]
#
## solver ------------------------------------------------------------------------------------
#        self.nl = self.conf.dict["Solver"]["Count"]
#

    def run(self):
        """
            run
            
            runs a simulation
            
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
#        t0, dt
#        y0 = [y01, ..., y0ny]
#        nl, nt, ny, nu, nb, nd, nexp, nimp
#        y, u, b, d, Aexp, Aimp
#        w
#        for l in range(nl):
#        for j in range(nt):
#            tj = t0 + j*dt
#
#            for ib in range(nb):
#                alpha, ialpha, beta, ibeta = _interpolate(nb[ib], tj)
#                bj[ib] = alpha * b[ib, ialpha] + beta * b[ib, ibeta]
#
#            for id in range(nd):
#                alpha, ialpha, beta, ibeta = _interpolate(nd[id], tj)
#                dj[id] = alpha * d[id, ialpha] + beta * d[id, ibeta]
#
#            alpha, ialpha, beta, ibeta = _interpolate(nexp, tj)
#            Aexpj = alpha * Aexp[ialpha] + beta * Aexp[ibeta]
#
#            alpha, ialpha, beta, ibeta = _interpolate(nimp, tj)
#            Aimpj = alpha * Aimpj[ialpha] + beta * Aimpj[ibeta]
#
#            qj = bgc(dt, tj, yj, u, bj, dj)
#            for i in range(ny):
#                Aexpj.multAdd(yj[i], qj[i], w)
#                Aimpj.mult(w, yj[i+1])
#
#            norm2 = norm2(yj[0] - yj[nt])


















