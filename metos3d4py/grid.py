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

from metos3d4py import util

# ----------------------------------------------------------------------------------------
class Grid:
    """
        grid context,
        
        grid masks, 2d/3d,
        compute vector length, number of profiles, profile depths,
        permutation index from nc/3d to tmm/1d

        """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        text = ""
#        text = text + "Grid:\n"
#        text = text + "shape: {}\n".format(self.mask3d.shape)
#        text = text + "profile count: {}\n".format(self.np)
#        text = text + "vector length: {}".format(self.nv)
        text = text + "        shape: {}\n".format(self.mask3d.shape)
        text = text + "profile count: {}\n".format(self.np)
        text = text + "vector length: {}".format(self.nv)
#        text = ""
#        text = text +
#        print("{}.__str__".format(self.__class__))
#        util.debug(m3d, self, "Shape: {}, profile count: {}, vector length: {}".format(mask3d.shape, np, nv), level=1)
#        return "Grid:\n  Vector length:   {}\n  Profile count:   {}\n  Profile lengths: {}".format(self.nv, self.np, self.npi)
        return text

# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):
        
        # get 'grid_mask' variable from file
        filepath = util.get_key(m3d, self, m3d.config, "Grid", str)
        gridfile = util.get_hdf5_file(m3d, self, filepath)
        grid = gridfile["grid_mask"]

        # masks
        self.mask3d = mask3d = (grid[...] != grid.fillvalue)
        self.mask2d = mask2d = mask3d[0,...]
        
        # not needed any more
        gridfile.close()

        # vector and profiles
        self.nv = nv = mask3d.sum()                     # vector length
        self.np = mask2d.sum()                          # profile count
        self.npi = mask3d.sum(axis=0)[mask2d]           # each profile length

        # index permutation from nc/3d to tmm/1d
        nc3d = mask3d[...].astype(">i4")
        (nz, ny, nx) = nc3d.shape
        nc3d[mask3d] = range(nv)
        nc3d[mask3d] = nc3d[mask3d] + 1
        nc1d = nc3d.reshape(nz, ny*nx).T.flat
        self.nc2tmm = nc1d[nc1d != 0] - 1

        # debug
        util.debug(m3d, self, self, level=1)

#        # store in self
#        self.mask3d = mask3d
#        self.mask2d = mask2d
#        self.nv = nv
#        self.np = np
#        self.npi = npi
#        self.nc2tmm = nc2tmm

#        util.debug(m3d, self, "Shape: {}, profile count: {}, vector length: {}".format(mask3d.shape, np, nv), level=1)

#        # store self in m3d
#        m3d.grid = self





