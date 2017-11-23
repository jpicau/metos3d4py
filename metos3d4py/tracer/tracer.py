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
        
        comm = m3d.comm
        grid = m3d.grid
        load = m3d.load
        
        # check tracer key
        try:
            conf = m3d.conf["Tracer"]
        except Exception as e:
            util._print_error(comm, "Tracer: Cannot retrieve tracer key from configuration.")
            sys.exit(1)

        tracer = conf["Name, Value, Unit, Description"]
        output = conf["Output"]

        input = conf.get("Input")
        if input is not None:
            util._print(comm, "Tracer init: Using input file: {}".format(input))
        
            try:
                file = h5py.File(input, "r")
            except Exception as e:
                util._print_error(comm, "Cannot open file: {}".format(input))
                util._print_error(comm, e)
                sys.exit(1)

        else:
            values = [t[1] for t in tracer]
            util._print(comm, "Tracer init: Using init values: {}".format(values))
                
        nv = grid.nv
        nvloc = load.nvloc
        
        ny = len(tracer)
        y0 = []
        for i in range(ny):
            y = PETSc.Vec()
            y.create()
            y.setType(PETSc.Vec.Type.STANDARD)
            y.setSizes((nvloc, nv))
            if input is not None:
                # read in
                util.read_from_nc_file(m3d, y, file, tracer[i][0], None)  # tracer name
#                util.set_from_nc_file(comm, grid, y, input, tracer[i][0], None)  # tracer name
            else:
                # set to tracer value
                y.set(tracer[i][1])
            y.assemble()
            y0.append(y)

        self.ny = ny
        self.y0 = y0
        self.output = output
        
        m3d.tracer = self

#        self.tracer = tracer
#        self.input = input
#        self.output = output
#        self.ny = ny
#        self.y = y
#        self.yj = yj
#        self.qj = qj
#        self.w = w





