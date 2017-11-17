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
    
    """

import sys
import yaml
import math
from ._version import version
from petsc4py.PETSc import COMM_WORLD as comm
from petsc4py import PETSc

# ----------------------------------------------------------------------------------------

def _print_usage():
    if comm.rank == 0:
        print("usage:\n  python metos3d.py [conf-yaml-file]")
        print("example:\n  python metos3d.py test/test.mitgcm-128x64x15.conf.yaml")

# ----------------------------------------------------------------------------------------

def _print_error(msg):
    if comm.rank == 0:
        print("### ERROR ### {}".format(msg))

# ----------------------------------------------------------------------------------------

def _print_message(msg):
    if comm.rank == 0:
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
        if comm.size > 1:
            _print_message("parallel run, {} processes".format(comm.size))
        else:
            _print_message("sequential run, {} process".format(comm.size))

        self.conf   = _parse_conf_yaml_file(argv)
        _print_message(self.conf)

# grid


        self.nl = self.conf["Solver"]["Count"]

#        t0, nt, dt = _get_time_step_conf()
        self.t0 = self.conf["Time step"]["Start"]
        self.nt = self.conf["Time step"]["Count"]
        self.dt = self.conf["Time step"]["Step"]
            
###        y0, ny = _get_tracer_conf()
##        self.y0 = y0
##        self.ny = ny
##        self.y0, ny
##          self.yj, qj, w
#        w = y0[0].duplicate()
#        for i in rnage(ny):
#            yj[i] = y0[i].duplicate()
#            qj[i] = y0[i].duplicate()
#
##self.Aexp, nexp
##    self.Aexpj
#        Aexp, nexp = _get_transport_explicit_conf()
#        Aexpj = Aexp.duplicate()
#
##self.Aimp, nimp
##    self.Aimpj
#        Aimp, nimp = _get_transport_implicit_conf()
#        Aimpj = Aimp.duplicate()
#
##self.u, nu
#        u, nu = _get_parameter_conf()
#
##self.b, nb, nbi
##    self.bj
#        b, nb, nbi = _get_boundary_data_conf()
#        for i in range(nb):
#            bj[i] = b[i][0].duplicate()
#
##self.d, nd, ndi
##    self.dj
#        d, nd, ndi = _get_domain_data_conf()
#        for i in range(nd):
#            dj[i] = d[i][0].duplicate()

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















#    def __init__(self):
#        """
#            __init__, constructor
#
#            initilize attributes to None,
#            print version and process count,
#
#        """
#        self.conf = None
#
#        _print_message("metos3d version {} ".format(version))
#        if comm.size > 1:
#            _print_message("parallel run, {} processes".format(comm.size))
#        else:
#            _print_message("sequential run, {} process".format(comm.size))
