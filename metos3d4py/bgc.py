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
from metos3d4py import util

class BGC:
    """
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
        return "BGC:\n  {}".format("...")
 
 # ----------------------------------------------------------------------------------------
    def init(self, m3d):
#        util.debug(m3d, "BGC init: {}".format("..."), level=1)
        util.debug(m3d, self, "{}".format("..."), level=1)

        try:
            config = m3d.config["BGC"]
        except Exception as e:
            util.debug(m3d, self, "No 'BGC' key found. Procceding without BGC model.", level=1)
            m3d.bgc = None
            return

#            # set to empty
#            self.nu = 0
#            self.u = []
#
#            self.nb = 0
#            self.nbi = []
#            self.b = []
##            self.bj = bj
#
#            self.nd = 0
#            self.ndi = []
#            self.d = []
##            self.dj = dj

#        self.init_parameter(m3d)
#        self.init_boundary_data(m3d)
#        self.init_domain_data(m3d)
        m3d.bgc = self

# ----------------------------------------------------------------------------------------
    def init_parameter(self, m3d):
        
        conf_bgc = m3d.config["BGC"]

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

#        self.domain_list = domain_list
#        self.domain_path = domain_path
    
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
    




