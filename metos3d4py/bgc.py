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

            y
            yj
            qj
            w
        
    """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        
        if not self.use_bgc:
            return "none"
                
        text = ""
        text = text + "name: {}\n".format(self.name)

        text = text + "parameter:\n"
        if self.use_param:
            text = text + "  {:<3} {:16} {:16} {:16} {:16}\n".format("no", "name", "value", "unit", "description")
            parameter = self.parameter
            nu = self.nu
            for i in range(nu):
                name = str(parameter[i][0])
                value = float(parameter[i][1])
                unit = str(parameter[i][2])
                descr = str(parameter[i][3])
                text = text + "  {:<3d} {:16.16} {:<16e} {:16.16} {:42.42}\n".format(i, name, value, unit, descr)
        else:
            text = text + "  none\n"

        text = text + "boundary data:\n"
        if self.use_boundary:
            text = text + "path: {}\n".format(self.boundary_path)
            text = text + "  {:<3} {:16} {:5} {:24} {:16} {:16}\n".format("no", "name", "count", "description", "unit", "file")
            boundary = self.boundary
            nb = self.nb
            for i in range(nb):
                name = str(boundary[i][0])
                count = self.nbi[i]
                descr = str(boundary[i][2])
                unit = str(boundary[i][3])
                file = str(boundary[i][4])
                text = text + "  {:<3d} {:16.16} {:>5d} {:24.24} {:16.16} {:16.16}\n".format(i, name, count, descr, unit, file)
        else:
            text = text + "  none\n"

        text = text + "domain data:\n"
        if self.use_domain:
            text = text + "path: {}\n".format(self.domain_path)
            text = text + "  {:<3} {:16} {:5} {:24} {:16} {:16}\n".format("no", "name", "count", "description", "unit", "file")
            domain = self.domain
            nd = self.nd
            for i in range(nd):
                name = str(domain[i][0])
                count = self.nbi[i]
                descr = str(domain[i][2])
                unit = str(domain[i][3])
                file = str(domain[i][4])
                text = text + "  {:<3d} {:16.16} {:>5d} {:24.24} {:16.16} {:16.16}\n".format(i, name, count, descr, unit, file)
        else:
            text = text + "  none\n"

        text = text.rstrip()

        return text

 # ----------------------------------------------------------------------------------------
    def __init__(self, m3d):
        
        self.config = config = m3d.config.get("BGC")
        if config is None:
            self.use_bgc = False
        else:
            self.use_bgc = True
            self.name = util.get_key(m3d, self, config, "Name", str)
            
            self.init_parameter(m3d)
            self.init_boundary_data(m3d)
            self.init_domain_data(m3d)

        # debug
        util.debug(m3d, self, self, level=1)

# ----------------------------------------------------------------------------------------
    def init_parameter(self, m3d):

        config = self.config.get("Parameter")
        if config is None:
            self.use_param = False
            return
        else:
            self.use_param = True
            
        self.parameter = parameter = util.get_key(m3d, self, config, "Name, Value, Unit, Description", list)

        self.nu = nu = len(parameter)
        u = []
        for i in range(nu):
            value = parameter[i][1]
            u.append(value)
        self.u = numpy.array(u)

# ----------------------------------------------------------------------------------------
    def init_boundary_data(self, m3d):
        
        config = self.config.get("Boundary data")
        if config is None:
            self.use_boundary = False
            return
        else:
            self.use_boundary = True
            
        self.boundary = boundary = util.get_key(m3d, self, config, "Name, Count, Description, Unit, File", list)
        self.boundary_path = util.get_key(m3d, self, config, "Path", str)

        np = m3d.grid.np
        nploc = m3d.load.nploc

        self.nb = nb = len(boundary)
        nbi = []
        b = []
        for i in range(nb):
            
            count = boundary[i][1]
            nbi.append(count)

            bi = util.create_petsc_vectors(nbi[i], (nploc, np))
            b.append(bi)

        self.nbi = nbi
        self.b = b
            
# ----------------------------------------------------------------------------------------
    def init_domain_data(self, m3d):

        config = self.config.get("Domain data")
        if config is None:
            self.use_domain = False
            return
        else:
            self.use_domain = True

        self.domain = domain = util.get_key(m3d, self, config, "Name, Count, Description, Unit, File", list)
        self.domain_path = util.get_key(m3d, self, config, "Path", str)

        nv = m3d.grid.nv
        nvloc = m3d.load.nvloc
        
        self.nd = nd = len(domain)
        ndi = []
        d = []
        for i in range(nd):
            
            count = domain[i][1]
            ndi.append(count)
            
            di = util.create_petsc_vectors(ndi[i], (nvloc, nv))
            d.append(di)

        self.ndi = ndi
        self.d = d

# ----------------------------------------------------------------------------------------
    def set(self):
        util.debug(m3d, self, "BGC set ...", level=1)

# ----------------------------------------------------------------------------------------
    def bgc(self, m3d, dt, qj, tj, yj, u, bj, dj):
        comm = m3d.comm
        conf_bgc = m3d.conf["BGC"]
        
        name = conf_bgc["Name"]
        util._print(comm, "BGC model: {}".format(name))
    











#            filepath = boundary_path + boundary_list[i][4]                              # file
#            util._print(comm, "BGC init boundary: Using input file: {}".format(filepath))
#
#            try:
#                file = h5py.File(filepath, "r")
#            except Exception as e:
#                util._print_error(comm, "Cannot open file: {}".format(filepath))
#                util._print_error(comm, e)
#                sys.exit(1)


#            for j in range(nbi[i]):
#                bw = PETSc.Vec()
#                bw.create()
#                bw.setType(PETSc.Vec.Type.STANDARD)
#                bw.setSizes((nploc, np))
#                # !!! read in !!!
#                util.read_from_nc_file(m3d, bw, file, boundary_list[i][0], j)     # boundary data name
#                # !!! read in !!!
#                bw.assemble()
#                b[i].append(bw)
#            # vector for interpolation
#            bj.append(b[i][0].duplicate())
#
#            file.close()

#            filepath = domain_path + domain_list[i][4]                                  # file
#            util._print(comm, "BGC init domain: Using input file: {}".format(filepath))
#
#            try:
#                file = h5py.File(filepath, "r")
#            except Exception as e:
#                util._print_error(comm, "Cannot open file: {}".format(filepath))
#                util._print_error(comm, e)
#                sys.exit(1)
#
#            for j in range(ndi[i]):
#                dw = PETSc.Vec()
#                dw.create()
#                dw.setType(PETSc.Vec.Type.STANDARD)
#                dw.setSizes((nvloc, nv))
#                # !!! read in !!!
#                util.read_from_nc_file(m3d, dw, file, domain_list[i][0], j)             # domain data name
#                # !!! read in !!!
#                dw.assemble()
#                d[i].append(dw)
#            dj.append(d[i][0].duplicate())
#
#            file.close()
#
##        self.domain_list = domain_list
##        self.domain_path = domain_path
#
#        self.nd = nd
#        self.ndi = ndi
#        self.d = d
#        self.dj = dj
