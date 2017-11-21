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

import h5py

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
    def init(self, comm, conf):
        file_path = conf.dict["Grid file"]
        grid_file = h5py.File(file_path, "r")
        grid = grid_file["grid_mask"]
        
        mask3d = (grid[...] != grid.fillvalue)
        mask2d = mask3d[0,...]
        
        nv = mask3d.sum()                           # vector length
        np = mask2d.sum()                           # profile count
        npi = mask3d.sum(axis=0)[mask2d]            # each profile length
        
        # compute index permutation from nc/3d to tmm/1d
        nc3d = mask3d[...].astype(">i4")
        (nz, ny, nx) = nc3d.shape
        nc3d[mask3d] = range(nv)
        nc3d[mask3d] = nc3d[mask3d] + 1
        nc1d = nc3d.reshape(nz, ny*nx).T.flat
        self.nc2tmm = nc1d[nc1d != 0] - 1

#        # compute index from nc/2d to tmm/boundary
#        nc2d = mask3d[0,...].astype(">i4")
#        nc2d[mask2d] = range(np)
#        nc2d[mask2d] = nc2d[mask2d] + 1

#        # test
#        test3d = mask3d[...].astype(">f8")
#        test3d[...] = 0.0
#        for iy in range(ny):
#            for ix in range(nx):
#                test3d[0, iy, ix] = float(iy) + float(ix)/1000.0
#        self.test3d = test3d

        self.mask3d = mask3d
        self.mask2d = mask2d
        self.nv = nv
        self.np = np
        self.npi = npi

    def __str__(self):
        return "Grid:\n  Vector length:   {}\n  Profile count:   {}\n  Profile lengths: {}".format(self.nv, self.np, self.npi)


