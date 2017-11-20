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

import numpy
from petsc4py import PETSc

class Bgc:
    """
        Bgc class
        
        biogeochemical model context
        
        Attributes:
            name
            tracer
            input
            output
            
            ny
            
            y
            yj
            qj
            w
        
        """

# ----------------------------------------------------------------------------------------

    def _init_tracer(self, conf, grid, load):
        
        name = conf.dict["Bgc model"]["Name"]
        tracer = conf.dict["Bgc model"]["Tracer"]["Name, Value, Unit, Description"]
        input = conf.dict["Bgc model"]["Tracer"]["Input"]
        output = conf.dict["Bgc model"]["Tracer"]["Output"]
        
        nv = grid.nv
        nvloc = load.nvloc
        
        ny = len(tracer)
        y = []
        yj = []
        qj = []
        for i in range(ny):
            yw = PETSc.Vec()
            yw.create()
            yw.setType(PETSc.Vec.Type.STANDARD)
            yw.setSizes((nvloc, nv))
            yw.assemble()
            # !!! read in !!!
            # !!! read in !!!
            y.append(yw)
            yj.append(yw.duplicate())
            qj.append(yw.duplicate())
        w = y[0].duplicate()
        
        self.name = name
        self.tracer = tracer
        self.input = input
        self.output = output
        
        self.ny = ny
        self.y = y
        self.yj = yj
        self.qj = qj
        self.w = w
    
# ----------------------------------------------------------------------------------------

    def _init_parameter(self, conf):

        # list of parameters
        parameter = conf.dict["Bgc model"]["Parameter"]["Name, Value, Unit, Description"]
        
        nu = len(parameter)
        u = []
        for i in range(nu):
            u.append(parameter[i][1])   # value
        u = numpy.array(u)

        self.parameter = parameter
        self.nu = nu
        self.u = u

# ----------------------------------------------------------------------------------------

    def _init_boundary_data(self, conf, grid, load):
        
        # list of boundary data
        boundary_list = conf.dict["Bgc model"]["Boundary data"]["Name, Count, Description, Unit, File"]
        boundary_path = conf.dict["Bgc model"]["Boundary data"]["Path"]

        np = grid.np
        nploc = load.nploc

        nb = len(boundary_list)
        b = []
        bj = []
        nbi = []
        for i in range(nb):
            b.append([])
            nbi.append(boundary_list[i][1])     # count
            for j in range(nbi[i]):
                bw = PETSc.Vec()
                bw.create()
                bw.setType(PETSc.Vec.Type.STANDARD)
                bw.setSizes((nploc, np))
                bw.assemble()
                # !!! read in !!!
                # !!! read in !!!
                b[i].append(bw)
            bj.append(b[i][0].duplicate())

        self.boundary_list = boundary_list
        self.boundary_path = boundary_path
        
        self.nb = nb
        self.nbi = nbi
        self.b = b
        self.bj = bj
            
# ----------------------------------------------------------------------------------------
    def _init_domain_data(self, conf, grid, load):
        
        # list of domain data
        domain_list = conf.dict["Bgc model"]["Domain data"]["Name, Count, Description, Unit, File"]
        domain_path = conf.dict["Bgc model"]["Domain data"]["Path"]
        
        nv = grid.nv
        nvloc = load.nvloc
        
        nd = len(domain_list)
        d = []
        dj = []
        ndi = []
        
        for i in range(nd):
            d.append([])
            ndi.append(domain_list[i][1])       # count
            for j in range(ndi[i]):
                dw = PETSc.Vec()
                dw.create()
                dw.setType(PETSc.Vec.Type.STANDARD)
                dw.setSizes((nvloc, nv))
                dw.assemble()
                # !!! read in !!!
                # !!! read in !!!
                d[i].append(dw)
            dj.append(d[i][0].duplicate())

        self.domain_list = domain_list
        self.domain_path = domain_path
    
        self.nd = nd
        self.ndi = ndi
        self.d = d
        self.dj = dj

# ----------------------------------------------------------------------------------------

    def init(self, comm, conf, grid, load):
        
        self._init_tracer(conf, grid, load)
        self._init_parameter(conf)
        self._init_boundary_data(conf, grid, load)
        self._init_domain_data(conf, grid, load)

# ----------------------------------------------------------------------------------------

    def __str__(self):
        return "Bgc:\n  Tracer name:  {}\n  Tracer count: {}".format(self.name, self.ny)
        


