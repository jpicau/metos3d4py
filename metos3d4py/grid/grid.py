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

import sys
import h5py
from metos3d4py.util import util

"""
    Grid submodule
    ==============
    
    Attributes:
        mask3d
        mask2d
        nv
        np
        npi
    
    """
class Grid:
    """
        an intance of the grid class is used to store the grid context,
        
        grid masks, 2d/3d,
        compute vector length, number of profiles, profile depths,
        permutation index from nc/3d to tmm/1d

        """

    def __str__(self):
        return "Grid:\n  Vector length:   {}\n  Profile count:   {}\n  Profile lengths: {}".format(self.nv, self.np, self.npi)

    def init(self, m3d):
        
        comm = m3d.comm
        conf = m3d.conf

        # check grid key
        try:
            file_path = conf["Grid"]
            util._print(comm, "Grid init: Using input file: {}".format(file_path))
        except Exception as e:
            util._print_error(comm, "Grid init: Cannot retrieve grid key from configuration.")
            sys.exit(1)
        
        # grid file
        try:
            grid_file = h5py.File(file_path, "r")
        except Exception as e:
            util._print_error(comm, "Grid init: Cannot open file.")
            util._print_error(comm, e)
            sys.exit(1)
        
        grid = grid_file["grid_mask"]
        
        # masks
        mask3d = (grid[...] != grid.fillvalue)
        mask2d = mask3d[0,...]
        
        # vector and profiles
        nv = mask3d.sum()                           # vector length
        np = mask2d.sum()                           # profile count
        npi = mask3d.sum(axis=0)[mask2d]            # each profile length
        
        # index permutation from nc/3d to tmm/1d
        nc3d = mask3d[...].astype(">i4")
        (nz, ny, nx) = nc3d.shape
        nc3d[mask3d] = range(nv)
        nc3d[mask3d] = nc3d[mask3d] + 1
        nc1d = nc3d.reshape(nz, ny*nx).T.flat
        nc2tmm = nc1d[nc1d != 0] - 1

        # store in self
        self.mask3d = mask3d
        self.mask2d = mask2d
        self.nv = nv
        self.np = np
        self.npi = npi
        self.nc2tmm = nc2tmm
    
        # store self in m3d
        m3d.grid = self



