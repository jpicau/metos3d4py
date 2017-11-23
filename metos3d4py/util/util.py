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

# ----------------------------------------------------------------------------------------
def usage(m3d):
    if m3d.rank == 0:
        print("usage:\n  python metos3d.py [conf-yaml-file]")
        print("example:\n  python metos3d.py test/test.mitgcm-128x64x15.conf.yaml")
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def error(m3d, msg):
    if m3d.rank == 0:
        for msg_line in str(msg).split("\n"):
            print("### ERROR ### {}".format(msg_line))
            sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def debug(m3d, msg, level=0):
    if m3d.rank == 0:
        if m3d.debug >= level:
            print(msg)
            sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def debug_synch(m3d, msg, level=0):
    """
        Print messages from each process in an ordered/synchronized manner.
        Collect each message first. Use the internal mpi4py communicator therefore.
        Rank 0 receives, all other process send.
    """
    comm = m3d.comm.tompi4py()
    size = m3d.size
    rank = m3d.rank
    
    if size > 1:
        if rank == 0:
            msgs = []
            msgs.append(str(msg))
            for i in range(size-1):
                req = comm.irecv(source=i+1, tag=i+1)
                msgs.append(str(req.wait()))
            debug(m3d, "\n".join(msgs), level=level)
        else:
            req = comm.isend(msg, dest=0, tag=rank)
            req.wait()
    else:
        debug(m3d, msg, level=level)

# ----------------------------------------------------------------------------------------
def get_key(m3d, caller, dict, key, valuetype):
    try:
        value = dict[key]
        if isinstance(value, valuetype):
            return value
        else:
            error(
                  m3d,
                  "{}: Value of key '{}' has type '{}'. Should be: '{}'"
                  .format(caller, key, type(value).__name__, valuetype.__name__))
            usage(m3d)
            sys.exit(1)
    except Exception as e:
        error(m3d, "{}: Cannot retrieve '{}' key from configuration.".format(caller, key))
        usage(m3d)
        sys.exit(1)

# ----------------------------------------------------------------------------------------
def get_file(m3d, caller, filepath):
    try:
        return open(filepath, "r")
    except Exception as e:
        error(m3d, "{}: Cannot open file: {}".format(caller, filepath))
        error(m3d, e)
        sys.exit(1)

# ----------------------------------------------------------------------------------------
def get_hdf5_file(m3d, caller, filepath):
    try:
        return h5py.File(filepath, "r")
    except Exception as e:
        error(m3d, "{}: Cannot open HDF5 file: {}".format(caller, filepath))
        error(m3d, e)
        sys.exit(1)

# ----------------------------------------------------------------------------------------
def set_vector_from_hdf5_file(m3d, v, file, varname, index):
    
    comm = m3d.comm
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
            _print_error(comm, "Variable: '{}' required to be 2D or 3D. Shape is: {} With index: {}".format(varname, var.shape, index))
            sys.exit(1)
    else:
        if len(var.shape) == 2:
            v[start:end] = var[...][mask2d][start:end]
        elif len(var.shape) == 3:
            v[start:end] = var[...][mask3d][nc2tmm][start:end]
        else:
            _print_error(comm, "Variable: '{}' required to be 2D or 3D. Shape is: {} ".format(varname, var.shape))
            sys.exit(1)

# ----------------------------------------------------------------------------------------
def get_config_from_yaml_file(m3d, argv):
    """
        parameters:
            argv        # command line arguments
        
        the first command line argument is always the name of the current executable,
        expect the file path as second argument,

        return:
            config      # configuration, dictionary where the parsed YAML contents is stored
    """

    if len(argv) > 1:
        f = get_file(m3d, "get_config_from_yaml_file", argv[1])
        
        try:
            # parse yaml content
            config = yaml.load(f)
        except Exception as e:
            error(m3d, "Cannot parse as YAML file: {}".format(argv[1]))
            error(m3d, e)
            sys.exit(1)

        f.close()
        return config

    else:
        error(m3d, "No configuration file given.")
        usage(m3d)
        sys.exit(0)

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






