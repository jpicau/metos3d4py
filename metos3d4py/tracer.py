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
from metos3d4py import util

class Tracer:
    """
        Tracer context
        
        Attributes:
            ny      number of tracers
            y0      initial tracers

        """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        text = "{:<3} {:16} {:16} {:16} {:16}\n".format("no", "name", "value", "unit","description")
        tracer = self.tracer
        ny = self.ny
        for i in range(ny):
            name = str(tracer[i][0])
            value = float(tracer[i][1])
            unit = str(tracer[i][2])
            descr = str(tracer[i][3])
            text = text + "{:<3d} {:16.16} {:<16e} {:16.16} {:42.42}\n".format(i, name, value, unit, descr)
        text = text + " input file: {}\n".format(self.input)
        text = text + "output file: {}".format(self.output)
        return text

# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):
        
        config = util.get_key(m3d, self, m3d.config, "Tracer", dict)
        self.tracer = tracer = util.get_key(m3d, self, config, "Name, Value, Unit, Description", list)
        self.output = config.get("Output")
        self.input = config.get("Input")

        # global and local vector length
        nv = m3d.grid.nv
        nvloc = m3d.load.nvloc
        
        self.ny = ny = len(tracer)
        self.y0 = util.create_petsc_vectors(ny, (nvloc, nv))

        # debug
        util.debug(m3d, self, self, level=1)

# ----------------------------------------------------------------------------------------
    def set(self, m3d):

        ny = self.ny
        y0 = self.y0
        tracer = self.tracer
        input = self.input

        text = ""
        if input is not None:
            file = util.get_hdf5_file(m3d, self, input)
            text = text + "set from input file:\n  {}\n".format(file.filename)
        else:
            values = [t[1] for t in tracer]
            text = text + "set from init values:\n"
            text = text + "  {:10} {}\n".format("name", "value")

        for i in range(ny):
            # Name, Value, Unit, Description
            name = str(tracer[i][0])
            value = float(tracer[i][1])
            if input is not None:
                util.set_vector_from_hdf5_file(m3d, y0[i], file, name)
            else:
                y0[i].set(value)
                text = text + "  {:10.10} {:<16e}\n".format(name, value)

        text = text.rstrip()
        util.debug(m3d, self, text, level=1)

# ----------------------------------------------------------------------------------------
    def save(self, m3d):
        util.debug(m3d, self, "Tracer: save {}".format(""), level=1)



