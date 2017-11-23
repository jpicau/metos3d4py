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
from petsc4py import PETSc
from metos3d4py.util import util

class Tracer:
    """
        Tracer context
        
        Attributes:
            y0      initial tracers

        """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        return "Tracer: {}".format("...")

# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        
        util.debug(m3d, "Tracer init: {}".format("..."), level=1)

        config = util.get_key(m3d, "Tracer init", m3d.config, "Tracer", dict)
        tracer = util.get_key(m3d, "Tracer init", config, "Name, Value, Unit, Description", list)
        output = util.get_key(m3d, "Tracer init", config, "Output", str)

        input = config.get("Input")
        if input is not None:
            util.debug(m3d, "Tracer init: Using input file: {}".format(input), level=1)
            tracerfile = util.get_hdf5_file(m3d, "Tracer init", input)
        else:
            tracervalue = [t[1] for t in tracer]
            util.debug(m3d, "Tracer init: Using init values: {}".format(tracervalue), level=1)

        # global and local vector length
        nv = m3d.grid.nv
        nvloc = m3d.load.nvloc
                
        # create tracer vectors
        ny = len(tracer)
        y0 = []
        for i in range(ny):
            y = PETSc.Vec()
            y.create()
            y.setType(PETSc.Vec.Type.STANDARD)
            y.setSizes((nvloc, nv))
            if input is not None:
                tracername = tracer[i][0]
                util.set_vector_from_hdf5_file(m3d, y, tracerfile, tracername, None)
            else:
                y.set(tracervalue[i])
            y.assemble()
            y0.append(y)

        self.ny = ny
        self.y0 = y0
        self.output = output

        m3d.tracer = self





