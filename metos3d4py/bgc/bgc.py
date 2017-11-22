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
import h5py
from petsc4py import PETSc
from metos3d4py.util import util

class BGC:
    """
        BGC class
        
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
    def __str__(self):
        return "BGC:\n  Model name:   {}\n  Tracer count: {}".format(self.name, self.ny)

# ----------------------------------------------------------------------------------------
    def init_tracer(self, m3d):
        
        comm = m3d.comm
        grid = m3d.grid
        load = m3d.load
        conf_bgc = m3d.conf["BGC"]
        
        name = conf_bgc["Name"]
        tracer = conf_bgc["Tracer"]["Name, Value, Unit, Description"]
        output = conf_bgc["Tracer"]["Output"]

        input = conf_bgc["Tracer"].get("Input")
        if input is not None:
            util._print(comm, "BGC init tracer: Using input file: {}".format(input))
        else:
            values = [t[1] for t in tracer]
            util._print(comm, "BGC init tracer: Using init values: {}".format(values))
                
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
            if input is not None:
                # read in
                util.set_from_nc_file(comm, grid, yw, input, tracer[i][0], None)  # tracer name
            else:
                # set to tracer value
                yw.set(tracer[i][1])
            yw.assemble()
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
    def init_parameter(self, m3d):
        
        conf_bgc = m3d.conf["BGC"]

        # list of parameters
        parameter = conf_bgc["Parameter"]["Name, Value, Unit, Description"]
        
        nu = len(parameter)
        u = []
        for i in range(nu):
            u.append(parameter[i][1])   # value
        u = numpy.array(u)

        self.parameter = parameter
        self.nu = nu
        self.u = u

# ----------------------------------------------------------------------------------------
    def init_boundary_data(self, m3d):
        
        comm = m3d.comm
        grid = m3d.grid
        load = m3d.load
        conf_bgc = m3d.conf["BGC"]

        # list of boundary data
        boundary_list = conf_bgc["Boundary data"]["Name, Count, Description, Unit, File"]
        boundary_path = conf_bgc["Boundary data"]["Path"]

        np = grid.np
        nploc = load.nploc

        nb = len(boundary_list)
        b = []
        bj = []
        nbi = []
        for i in range(nb):
            b.append([])
            nbi.append(boundary_list[i][1])                                             # count

            filepath = boundary_path + boundary_list[i][4]                              # file
            util._print(comm, "BGC init boundary: Using input file: {}".format(filepath))

            try:
                file = h5py.File(filepath, "r")
            except Exception as e:
                util._print_error(comm, "Cannot open file: {}".format(filepath))
                util._print_error(comm, e)
                sys.exit(1)
            
            for j in range(nbi[i]):
                bw = PETSc.Vec()
                bw.create()
                bw.setType(PETSc.Vec.Type.STANDARD)
                bw.setSizes((nploc, np))
                # !!! read in !!!
                util.read_from_nc_file(m3d, bw, file, boundary_list[i][0], j)     # boundary data name
                # !!! read in !!!
                bw.assemble()
                b[i].append(bw)
            # vector for interpolation
            bj.append(b[i][0].duplicate())
                
            file.close()

        self.boundary_list = boundary_list
        self.boundary_path = boundary_path
        
        self.nb = nb
        self.nbi = nbi
        self.b = b
        self.bj = bj
            
# ----------------------------------------------------------------------------------------
    def init_domain_data(self, m3d):

        comm = m3d.comm
        grid = m3d.grid
        load = m3d.load
        conf_bgc = m3d.conf["BGC"]

        # list of domain data
        domain_list = conf_bgc["Domain data"]["Name, Count, Description, Unit, File"]
        domain_path = conf_bgc["Domain data"]["Path"]
#        domain_list = conf["Name, Count, Description, Unit, File"]
#        domain_path = conf["Path"]

        nv = grid.nv
        nvloc = load.nvloc
        
        nd = len(domain_list)
        d = []
        dj = []
        ndi = []
        
        for i in range(nd):
            d.append([])
            ndi.append(domain_list[i][1])                                               # count
            
            filepath = domain_path + domain_list[i][4]                                  # file
            util._print(comm, "BGC init domain: Using input file: {}".format(filepath))

            try:
                file = h5py.File(filepath, "r")
            except Exception as e:
                util._print_error(comm, "Cannot open file: {}".format(filepath))
                util._print_error(comm, e)
                sys.exit(1)

            for j in range(ndi[i]):
                dw = PETSc.Vec()
                dw.create()
                dw.setType(PETSc.Vec.Type.STANDARD)
                dw.setSizes((nvloc, nv))
                # !!! read in !!!
                util.read_from_nc_file(m3d, dw, file, domain_list[i][0], j)             # domain data name
                # !!! read in !!!
                dw.assemble()
                d[i].append(dw)
            dj.append(d[i][0].duplicate())
                
            file.close()

        self.domain_list = domain_list
        self.domain_path = domain_path
    
        self.nd = nd
        self.ndi = ndi
        self.d = d
        self.dj = dj

# ----------------------------------------------------------------------------------------
    def bgc(self, m3d, dt, qj, tj, yj, u, bj, dj):
        comm = m3d.comm
        conf_bgc = m3d.conf["BGC"]
        
        name = conf_bgc["Name"]
        util._print(comm, "BGC model: {}".format(name))
    
    
# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        self.init_tracer(m3d)
        self.init_parameter(m3d)
        self.init_boundary_data(m3d)
        self.init_domain_data(m3d)




