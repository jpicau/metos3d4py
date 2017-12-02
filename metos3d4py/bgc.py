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

import importlib
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
            
            # load bgc module
            # must be placed in metos3d/bgc/<<name>>/<<sources>>
            # must be compiled with: f2py -c <<sources>> -m <<name>>
            # call to bgc routine must apply to:
            # fortran: metos3dbgc(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d)
            # python: metos3dbgc(dt,q,t,y,u,b,d,[ny,nx,nu,nb,nd])
            # module name: exchange '-' to '_'
            modname = self.name.replace("-", "_")
            self.mod = importlib.import_module("bgc." + self.name + "." + modname)

            self.init_parameter(m3d)
            self.init_boundary_data(m3d)
            self.init_domain_data(m3d)

        # debug
        util.debug(m3d, self, self, level=1)

# ----------------------------------------------------------------------------------------
    def init_parameter(self, m3d):

        self.u = numpy.array([], dtype="f8")    # f2py, no copy as such
        self.nu = 0

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
        
        self.b = []
        self.nb = 0
        self.nbi = []

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

        self.d = []
        self.nd = 0
        self.ndi = []

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
    def set(self, m3d):

        text = ""
        if self.use_boundary:
            nb = self.nb
            nbi = self.nbi
            b = self.b
            boundary = self.boundary
            path = self.boundary_path
            text = text + "boundary:\n"
            text = text + "  {:10.10} {}\n".format("name", "file")
            for i in range(nb):
                # Name, Count, Description, Unit, File
                varname = boundary[i][0]
                filename = boundary[i][4]
                file = util.get_hdf5_file(m3d, self, path + filename)
                text = text + "  {:10.10} {}\n".format(varname, path + filename)
                for j in range(nbi[i]):
                    pass
                    util.set_vector_from_hdf5_file(m3d, b[i][j], file, varname, index=j)

        if self.use_domain:
            nd = self.nd
            ndi = self.ndi
            d = self.d
            domain = self.domain
            path = self.domain_path
            text = text + "domain:\n"
            text = text + "  {:10.10} {}\n".format("name", "file")
            for i in range(nd):
                # Name, Count, Description, Unit, File
                varname = domain[i][0]
                filename = domain[i][4]
                file = util.get_hdf5_file(m3d, self, path + filename)
                text = text + "  {:10.10} {}\n".format(varname, path + filename)
                for j in range(ndi[i]):
                    util.set_vector_from_hdf5_file(m3d, d[i][j], file, varname, index=j)
        
        text = text.rstrip()
        util.debug(m3d, self, text, level=1)

## ----------------------------------------------------------------------------------------
#    def bgc(self, m3d, dt, qj, tj, yj, u, bj, dj):
#        comm = m3d.comm
#        conf_bgc = m3d.conf["BGC"]
#
#        name = conf_bgc["Name"]
#        util._print(comm, "BGC model: {}".format(name))

## ----------------------------------------------------------------------------------------
#    def evaluate(self, m3d, dt, qj, tj, yj, u, bj, dj):
#        pass

# ----------------------------------------------------------------------------------------
    def q(self, m3d, dt, q, t, y, u, b, d):
    
        ny = m3d.tracer.ny
        npi = m3d.grid.npi
        npiloc = m3d.load.npiloc
        nvprev = m3d.load.nvprev
        nploc = m3d.load.nploc
        npprev = m3d.load.npprev
        indptr = m3d.load.indptr
        indptrloc = m3d.load.indptrloc

        din = numpy.array([], dtype="f8")
        bin = numpy.array([], dtype="f8")

        # stack lists, access as numpy array then
        ys = numpy.stack(y)
        qs = numpy.stack(q)
#        qs = numpy.stack([qi.array_w[...] for qi in q])
        if d: ds = numpy.stack(d)
        if b: bs = numpy.stack(b)

#        print(yin[indptrloc[0]:indptrloc[1]])

#        print(indptr)
#        print(indptrloc)
#        print(npi)
#        print(npiloc)

#        print(type(q[0][indptr[0]:indptr[1]]))
#        print(q[:])
#        print([type(vi[indptr[0]:indptr[1]]) for vi in y])
#            start, end = A.owner_range

#        y[:][start:end]
#        ip = 0
#        # local arrays
#        yin = numpy.stack([yi.array for yi in y]).T
##        yin = numpy.stack([yi.array for yi in y])
#        print(yin)
#        print(yin.flags["F_CONTIGUOUS"])

#        print(numpy.split(y[0].array, indptr[1:-1]))

#        print(yin[indptr[ip]:indptr[ip+1],:].flags["F_CONTIGUOUS"])

#        print(zip(*[yi.array for yi in y]))

#        print(id(yin[indptr[ip]:indptr[ip+1],:]))
#        print(id(yin[indptr[ip]:indptr[ip+1]]))

#        print(npi[npprev:npprev+nploc])
#        print(id(y[0][indptr[0]:indptr[1]]))
#        print(id(y[0][...]))
#        print(id(y[0].array))

#        import sys
#        sys.stdout.flush()
#        sys.exit(1)

#        for ip in range(npprev, npprev+nploc):
        for ip in range(nploc):

#            q_inout = numpy.zeros((npi[npprev+ip], ny), order="F")
#            q_inout = numpy.zeros((npiloc[ip], ny), order="F")
#            q_inout_c = numpy.array(q_inout, order="C")
#            print(q_inout_c.shape)
#            print(q_inout.shape)
#            sys.stdout.flush()

#            # fortran contiguous arrays
#            yin = numpy.stack([yi[indptr[npprev+ip]:indptr[npprev+ip+1]] for yi in y]).T
#            if d:
#                din = numpy.stack([di[indptr[npprev+ip]:indptr[npprev+ip+1]] for di in d]).T
#            if b:
#                bin = numpy.stack([bi[npprev+ip] for bi in b]).T

#            print(ip)
#            print(indptr[ip])
#            print(yin[indptr[ip]:indptr[ip+1],:])
#            print(yin[nvprev+indptr[npprev+ip]:nvprev+indptr[npprev+ip+1]])

#            print(indptrloc[ip])
#            print(indptrloc[ip+1])
#            print(yin[0][0])
#            print(yin[:,indptrloc[ip]:indptrloc[ip+1]].T)

            import sys

#            print(ys.shape)
#            print(ds.shape)
#            print(bs.shape)
#            yin = ys[:,indptrloc[ip]:indptrloc[ip+1]]
#            yin = ys[:,indptrloc[ip]:indptrloc[ip+1]].T
            yin = numpy.asfortranarray(ys[:,indptrloc[ip]:indptrloc[ip+1]].T)
            qin = numpy.asfortranarray(qs[:,indptrloc[ip]:indptrloc[ip+1]].T)
#            print(qs[:,indptrloc[ip]:indptrloc[ip+1]].T.flags["F_CONTIGUOUS"])
            if d: din = numpy.asfortranarray(ds[:,indptrloc[ip]:indptrloc[ip+1]].T)
            if b: bin = numpy.asfortranarray(bs[:,ip])
#            yin = ys[:,indptr[ip]:indptr[ip+1]]
#            print((indptr[ip], indptr[ip+1]))
#            print((indptrloc[ip], indptrloc[ip+1]))
#            print(yin)
#            print(din)
#            print(bin)
#            print(qin)


#            import sys
            sys.stdout.flush()
#            sys.exit(1)
            # call fortran model
            # python: metos3dbgc(dt,q,t,y,u,b,d,[ny,nx,nu,nb,nd])
            self.mod.metos3dbgc(
                                dt,
#                                yin[indptr[ip]:indptr[ip+1],:],
#                                q[0][indptr[ip]:indptr[ip+1]],
#                                q_inout,
                                qin,
                                t,
                                yin,
#                                y[0][indptr[npprev+ip]:indptr[npprev+ip+1]],
#                                y[0][nvprev+indptr[ip]:nvprev+indptr[ip+1]],
#                                yin[indptr[ip]:indptr[ip+1],:],
                                u,
#                                b[0][npprev+ip],
#                                numpy.array([], dtype="f8"),
                                bin,
#                                d[0][indptr[npprev+ip]:indptr[npprev+ip+1]],
#                                numpy.array([], dtype="f8"),
                                din,
                                ny,
                                npiloc[ip],
                                self.nu,
                                self.nb,
                                self.nd)
#            print(qin)
            # copy back result, global indices
            for i in range(ny):
#                q[i][indptrloc[ip]:indptrloc[ip+1]] = qin[:,i]
                q[i][indptr[npprev+ip]:indptr[npprev+ip+1]] = qin[:,i]
#            print(qs[:,indptrloc[ip]:indptrloc[ip+1]].shape)
#            print(qs[:,indptrloc[ip]:indptrloc[ip+1]])
#            print(qin.shape)
#            print(qin)
#            qs[:,indptrloc[ip]:indptrloc[ip+1]] = qin.T
#            print(qs[:,indptrloc[ip]:indptrloc[ip+1]])
#            print(q[0][indptrloc[ip]:indptrloc[ip+1]].shape)
#            print(q[0][indptrloc[ip]:indptrloc[ip+1]])
##            print(q[:].array.shape)
#            print("-")
##            qin = numpy.asfortranarray(qs[:,indptrloc[ip]:indptrloc[ip+1]].T)
#
##            print(q_inout)
            sys.stdout.flush()
##
#            sys.exit(1)

        for qi in q:
            qi.assemble()









