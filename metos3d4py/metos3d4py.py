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

import sys
import yaml
#import math
import numpy
import h5py
from ._version import version
#from petsc4py.PETSc import COMM_WORLD as comm
from petsc4py import PETSc

# ----------------------------------------------------------------------------------------

def _print_usage():
    if PETSc.COMM_WORLD.rank == 0:
        print("usage:\n  python metos3d.py [conf-yaml-file]")
        print("example:\n  python metos3d.py test/test.mitgcm-128x64x15.conf.yaml")

# ----------------------------------------------------------------------------------------

def _print_error(msg):
    if PETSc.COMM_WORLD.rank == 0:
        print("### ERROR ### {}".format(msg))

# ----------------------------------------------------------------------------------------

def _print_message(msg):
    if PETSc.COMM_WORLD.rank == 0:
        print(msg)

# ----------------------------------------------------------------------------------------

def _parse_conf_yaml_file(argv):

    # the first command line argument is always the name of the current executable,
    # expect the file path as second argument,
    if len(argv) > 1:
        
        _print_message("parsing configuration file: {}".format(argv[1]))
        
        try:
            # open file
            f = open(argv[1])
        except Exception as e:
            _print_error("Cannot open file: {}".format(argv[1]))
            _print_message(e)
            _print_usage()
            sys.exit(1)
    
        try:
            # parse yaml content
            conf = yaml.load(f)
        except Exception as e:
            _print_error("Cannot parse as YAML file: {}".format(argv[1]))
            _print_message(e)
            _print_usage()
            sys.exit(1)
        
        f.close()
        return conf

    else:
        _print_error("No configuration file given.")
        _print_usage()
        sys.exit(1)

# ----------------------------------------------------------------------------------------

class Metos3D:
    """
        Metos3D class
        
        serves as the Metos3D data type, i.e.
        stores the metos3d simulation context,
        
        conf:       configuration dictionary
        
        init():     initialize metos3d context from a yaml configuration file
        run():      perform a simulation
        final():    finalize the context
    
    """
    def __init__(self):
        """
            __init__, constructor

            initilize attributes,

        """
        self.version = version
        self.comm = PETSc.COMM_WORLD
        self.size = self.comm.size
        self.rank = self.comm.rank
        
#        o = PETSc.Options()
##        o.setValue("-help", True)
#        o.setValue("-malloc_dump", True)
##        o.view()

        self.conf = None
    

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
            
            argv:   list with command arguments
            
            initialize metos3d context,
            do on every process the same,
            
            print version and process count,
            
            load and parse yaml configuration file,
            store configuration dictionary in instance variabel,
            """
        
        _print_message("metos3d version {} ".format(version))
        if self.comm.size > 1:
            _print_message("parallel run, {} processes".format(self.comm.size))
        else:
            _print_message("sequential run, {} process".format(self.comm.size))

        self.conf   = _parse_conf_yaml_file(argv)
#        _print_message(self.conf)

## grid ------------------------------------------------------------------------------------
#grid = self.Grid().loadFromFile(file_path)

        filepath = self.conf["Grid file"]
        gridfile = h5py.File(filepath, "r")
        
        grid = gridfile["grid_mask"]
#        (nz, ny, nx) = grid.shape
#        print("{} {} {}".format(nz, ny, nx))

#        print(grid[:]grid.fillvalue)
#        print(type(grid[...]))
        mask = (grid[...] != grid.fillvalue)
        mask2d = mask[0,...]
        
        self.nv = mask.sum()                     # vector length
        self.np = mask2d.sum()                   # profile count
        self.npi = mask.sum(axis=0)[mask2d]      # profile length each

#        print(nv)
#        print(np)
#        print(npi)

#        gridfile.close()

## load ------------------------------------------------------------------------------------
#load = self.Load().

        starts = self.npi.cumsum()-self.npi
        weights = starts + 0.5*self.npi
#        print(starts)
#        print(weights)

        ranks = numpy.floor((weights/self.nv)*self.size)
#        print(ranks)

        rank_unique, rank_index, rank_count = numpy.unique(ranks, return_index = True, return_counts = True)
        rank_prev = rank_count.cumsum() - rank_count
#        print((rank_unique, rank_index, rank_count))

        self.nploc = rank_count[self.rank]
        self.npprev = rank_prev[self.rank]

        self.nvloc = self.npi[self.npprev:self.npprev+self.nploc].sum()
        self.nvprev = self.npi[0:self.npprev].sum()

#        print("{:2} profile: {:5d} {:5d} vector: {:5d} {:5d}".format(
#            self.rank, self.npprev, self.nploc, self.nvprev, self.nvloc))

# bgc ------------------------------------------------------------------------------------
## tracer -------------------------------------

#        self.bgc = Bgc()
#        self.bgc.name = self.conf["Bgc model"]["Name"]
        tracer_list = self.conf["Bgc model"]["Tracer"]["Name, Value, Unit, Description"]
        self.ny = len(tracer_list)
        y = []
        yj = []
        qj = []
        for i in range(self.ny):
#            print("{} {}".format(i, tracer_list[i]))
            yw = PETSc.Vec()
            yw.create()
            yw.setType(PETSc.Vec.Type.STANDARD)
            yw.setSizes((self.nvloc, self.nv))
            yw.assemble()
            # !!! read in !!!
            # !!! read in !!!
            y.append(yw)
            yj.append(yw.duplicate())
            qj.append(yw.duplicate())
#        print(y0)
        self.y = y
        self.yj = yj
        self.qj = qj
        self.w = y[0].duplicate()
        
        self.model_name = self.conf["Bgc model"]["Name"]
        self.input_file = self.conf["Bgc model"]["Tracer"]["Input"]
        self.output_file = self.conf["Bgc model"]["Tracer"]["Output"]

###        y0, ny = _get_tracer_conf()
##        self.y0 = y0
##        self.ny = ny
##        self.y0, ny
##          self.yj, qj, w
#        w = y0[0].duplicate()
#        for i in rnage(ny):
#            yj[i] = y0[i].duplicate()
#            qj[i] = y0[i].duplicate()

## parameter -------------------------------------
##self.u, nu
#        u, nu = _get_parameter_conf()
        u = []
        parameter_list = self.conf["Bgc model"]["Parameter"]["Name, Value, Unit, Description"]
        nu = len(parameter_list)
        for i in range(nu):
#            print(parameter_list[i])
            u.append(parameter_list[i][1]) # value
        u = numpy.array(u)
#        print(u)

## forcing, boundary -------------------------------------
##self.b, nb, nbi
##    self.bj
#        b, nb, nbi = _get_boundary_data_conf()
#        for i in range(nb):
#            bj[i] = b[i][0].duplicate()
        b = []
        bj = []
        boundary_list = self.conf["Bgc model"]["Boundary data"]["Name, Count, Description, Unit, File"]
        boundary_path = self.conf["Bgc model"]["Boundary data"]["Path"]
        nb = len(boundary_list)
        nbi = []
        for i in range(nb):
#            print(boundary_list[i])
            b.append([])
            nbi.append(boundary_list[i][1])
#            print(nbi)
            for j in range(nbi[i]):
                bw = PETSc.Vec()
                bw.create()
                bw.setType(PETSc.Vec.Type.STANDARD)
                bw.setSizes((self.nploc, self.np))
                bw.assemble()
                # !!! read in !!!
                # !!! read in !!!
                b[i].append(bw)
            bj.append(b[i][0].duplicate())
#        print(b)

## forcing, domain -------------------------------------
##self.d, nd, ndi
##    self.dj
#        d, nd, ndi = _get_domain_data_conf()
#        for i in range(nd):
#            dj[i] = d[i][0].duplicate()
        domain_list = self.conf["Bgc model"]["Domain data"]["Name, Count, Description, Unit, File"]
        domain_path = self.conf["Bgc model"]["Domain data"]["Path"]
        d = []
        dj = []
        nd = len(domain_list)
        ndi = []
        for i in range(nd):
#            print(domain_list[i])
            d.append([])
            ndi.append(domain_list[i][1])
#            print(ndi)
            for j in range(ndi[i]):
                dw = PETSc.Vec()
                dw.create()
                dw.setType(PETSc.Vec.Type.STANDARD)
                dw.setSizes((self.nvloc, self.nv))
                dw.assemble()
                # !!! read in !!!
                # !!! read in !!!
                d[i].append(dw)
            dj.append(d[i][0].duplicate())
#        print(d)

# transport ------------------------------------------------------------------------------------
##self.Aexp, nexp
##    self.Aexpj
#        Aexp, nexp = _get_transport_explicit_conf()
#        Aexpj = Aexp.duplicate()
##self.Aimp, nimp
##    self.Aimpj
#        Aimp, nimp = _get_transport_implicit_conf()
#        Aimpj = Aimp.duplicate()
        transport_list = self.conf["Transport"]["Name, Count, Type, File"]
        transport_patj = self.conf["Transport"]["Path"]
        

# time step
#        t0, nt, dt = _get_time_step_conf()
        self.t0 = self.conf["Time step"]["Start"]
        self.nt = self.conf["Time step"]["Count"]
        self.dt = self.conf["Time step"]["Step"]

# solver ------------------------------------------------------------------------------------
        self.nl = self.conf["Solver"]["Count"]

# ----------------------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------------------

def _interpolate(n, t):
    '''
        _interpolate
        
        n:  number of intervals [0,1] is devided into
        t:  point in time (in [0,1])
        
        compute the interpolation coefficients and indices on the fly
        
        '''
    
    w = t * n + 0.5
    beta = math.fmod(w, 1.0)
    alpha = (1.0 - beta)
    ibeta = int(math.fmod(math.floor(w), n))
    ialpha = int(math.fmod(math.floor(w) + n - 1, n))
    
    return alpha, ialpha, beta, ibeta

# ----------------------------------------------------------------------------------------















#
#        _print_message("metos3d version {} ".format(version))
#        if comm.size > 1:
#            _print_message("parallel run, {} processes".format(comm.size))
#        else:
#            _print_message("sequential run, {} process".format(comm.size))

#        print(tuple(mask.sum(axis=0)[4,:]))
#        print("{:d}".format(mask.sum(axis=0)[4,:]))
#        print("%d"*128 % tuple(mask.sum(axis=0)[4,:]))
#        print("%d"*128 % mask.sum(axis=0)[4,:])
#        line = map("%d".format, mask.sum(axis=0)[4,:])
#        print(line)

#        print("%d"*128 % mask.sum(axis=0)[4,:])
#        print(str(mask.sum(axis=0)[4,:]))
#        print("{}".format(mask.sum(axis=0)[4,:]))
#        print(mask.sum(axis=0).flatten().nonzero())

#        print(grid[mask])
#        print(mask.sum())
#        print(mask[0,...].sum())

#        mask = mask.reshape(nz, ny*nx)
#        mask = mask.T

#        print(mask.flatten())
#        print(mask.sum(axis=1))
#        print(mask.size)
#        print(sum(sum(sum(~mask[:,:,:]))))
#        print(sum(sum(~mask[0,:,:])))

#        print(gridh5["grid_mask"])
#        print(dir(gridh5))
#        for k,v in gridh5.items():
#            print("{}, {}".format(k, v))

#        viewerhdf5 = PETSc.ViewerHDF5()
#        viewerhdf5.create(gridfile, mode = "r", comm = comm)
#        print(viewerhdf5)

