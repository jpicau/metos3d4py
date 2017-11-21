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

import math
import sys
import h5py
from petsc4py import PETSc

# ----------------------------------------------------------------------------------------
def _print_usage(comm):
    if comm.rank == 0:
        print("usage:\n  python metos3d.py [conf-yaml-file]")
        print("example:\n  python metos3d.py test/test.mitgcm-128x64x15.conf.yaml")
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_error(comm, msg):
    if comm.rank == 0:
        print("### ERROR ### {}".format(msg))
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_message(comm, msg):
    if comm.rank == 0:
        print(msg)
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_message_synch(comm, msg):
    
    comm_mpi = comm.tompi4py()
    size = comm.size
    rank = comm.rank
    
    if size > 1:
        if rank == 0:
            msgs = []
            msgs.append(str(msg))
            for i in range(size-1):
                req = comm_mpi.irecv(source=i+1, tag=i+1)
                msgs.append(str(req.wait()))
            _print_message(comm, "\n".join(msgs))
        else:
            req = comm_mpi.isend(msg, dest=0, tag=rank)
            req.wait()
    else:
        _print_message(comm, msg)

# ----------------------------------------------------------------------------------------
def _set_from_nc_file(comm, grid, v, varfile, varname, index):
    
    try:
        file = h5py.File(varfile, "r")
    except Exception as e:
        _print_error(comm, "Cannot open file: {}".format(varfile))
        _print_message(comm, e)
        sys.exit(1)

    try:
        var = file[varname]
    except Exception as e:
        _print_error(comm, "Cannot retrieve variable: {} from file: {}".format(varname, varfile))
        _print_message(comm, e)
        sys.exit(1)

#    _print_message(comm, "{} {} {}".format(varname, index, var.shape))

    start, end = v.getOwnershipRange()
#    print("{} {}".format(start, end))

    if index is not None:
        # C order, slowest dim left
        if len(var.shape) == 3:
            # mask2d
            v[start:end] = var[index,...][grid.mask2d][start:end]
        elif len(var.shape) == 4:
            # mask3d
            v[start:end] = var[index,...][grid.mask3d][grid.nc2tmm][start:end]
        else:
            _print_error(comm, "Variable: {} required to be 2D or 3D. Shape is: {} With index: {}".format(varname, var.shape, index))
    else:
        if len(var.shape) == 2:
            v[start:end] = var[...][grid.mask2d][start:end]
        elif len(var.shape) == 3:
            v[start:end] = var[...][grid.mask3d][grid.nc2tmm][start:end]
        else:
            _print_error(comm, "Variable: {} required to be 2D or 3D. Shape is: {} ".format(varname, var.shape))

    file.close()

#    _print_message_synch(comm, "{:3} {:7} {:7} {}".format(comm.rank, start, end, v[start:end][:4]))


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

