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

import yaml
import math
import sys
import h5py
from petsc4py import PETSc

"""
    read_conf_from_yaml_file
    
    """

# ----------------------------------------------------------------------------------------
def _print(comm, msg):
    if comm.rank == 0:
        print(msg)
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_usage(comm):
    if comm.rank == 0:
        print("usage:\n  python metos3d.py [conf-yaml-file]")
        print("example:\n  python metos3d.py test/test.mitgcm-128x64x15.conf.yaml")
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_error(comm, msg):
    if comm.rank == 0:
        for msg_line in str(msg).split("\n"):
            print("### ERROR ### {}".format(msg_line))
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def _print_synch(comm, msg):
    """
        Print messages from each process in an ordered/synchronized manner.
        Collect each message first. Use the internal mpi4py communicator therefore.
        Rank 0 receives, all other process send.
        
        """
    
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
def read_from_nc_file(m3d, v, file, varname, index):
    
    grid = m3d.grid
    mask3d = grid.mask3d
    mask2d = grid.mask2d
    nc2tmm = grid.nc2tmm

    try:
        var = file[varname]
    except Exception as e:
        _print_error(comm, "Cannot retrieve variable: {}".format(varname))
        _print_error(comm, e)
        sys.exit(1)

    start, end = v.getOwnershipRange()

    if index is not None:
        # C order, slowest dim left
        if len(var.shape) == 3:
            # mask2d
            v[start:end] = var[index,...][mask2d][start:end]
        elif len(var.shape) == 4:
            # mask3d
            v[start:end] = var[index,...][mask3d][nc2tmm][start:end]
        else:
            util._print_error(comm, "Variable: {} required to be 2D or 3D. Shape is: {} With index: {}".format(varname, var.shape, index))
    else:
        if len(var.shape) == 2:
            v[start:end] = var[...][mask2d][start:end]
        elif len(var.shape) == 3:
            v[start:end] = var[...][mask3d][nc2tmm][start:end]
        else:
            util._print_error(comm, "Variable: {} required to be 2D or 3D. Shape is: {} ".format(varname, var.shape))

# ----------------------------------------------------------------------------------------
def interp(n, t):
    '''
        compute the interpolation coefficients and indices on the fly

        n:  number of intervals [0,1] is devided into
        t:  point in time (in [0,1])
        
        alpha, ialpha, beta, ibeta
        
        '''

    w = t * n + 0.5
    beta = math.fmod(w, 1.0)
    alpha = (1.0 - beta)
    ibeta = int(math.fmod(math.floor(w), n))
    ialpha = int(math.fmod(math.floor(w) + n - 1, n))
    
    return alpha, ialpha, beta, ibeta

# ----------------------------------------------------------------------------------------
def read_conf_from_yaml_file(m3d, argv):
    """
        argv    # command line arguments
        conf    # dictionary where the parsed YAML contents is stored
            
        """
    
    comm = m3d.comm
    
    # the first command line argument is always the name of the current executable,
    # expect the file path as second argument,
    if len(argv) > 1:
        
        _print(comm, "parsing configuration file: {}".format(argv[1]))

        try:
            # open file
            f = open(argv[1])
        except Exception as e:
            _print_error(comm, "Cannot open file: {}".format(argv[1]))
            _print_error(comm, e)
            _print_usage(comm)
            sys.exit(1)

        try:
            # parse yaml content
            m3d.conf = yaml.load(f)
        except Exception as e:
            _print_error(comm, "Cannot parse as YAML file: {}".format(argv[1]))
            _print_error(comm, e)
            _print_usage(comm)
            sys.exit(1)

        f.close()
    
    else:
        _print_error(comm, "No configuration file given.")
        _print_usage(comm)
        sys.exit(0)

