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
import numpy
from petsc4py import PETSc
from metos3d4py import util

class TMM:
    """
        TMM class
        
        Transport matrix context
        
        Attributes:
            nexp
            Aexp
            Aexpj
            nimp
            Aimp
            Aimpj

        """
    
# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):

        self.config = config = m3d.config.get("TMM")
        if config is None:
            self.use_tmm = False
        else:
            self.use_tmm = True

            self.path = util.get_key(m3d, self, config, "Path", str)
            self.init_explicit(m3d)
            self.init_implicit(m3d)

        # debug
        util.debug(m3d, self, self, level=1)

# ----------------------------------------------------------------------------------------
    def init_explicit(self, m3d):
    
        self.explicit = explicit = self.config.get("Explicit")
        if explicit is None:
            self.use_explicit = False
            return
        else:
            self.use_explicit = True
        
# ----------------------------------------------------------------------------------------
    def init_implicit(self, m3d):

        self.implicit = implicit = self.config.get("Implicit")
        if implicit is None:
            self.use_implicit = False
            return
        else:
            self.use_implicit = True

# ----------------------------------------------------------------------------------------
    def set(self, m3d):
        
        if not self.use_tmm:
            util.debug(m3d, self, "none", level=1)
            return
        else:
            nv = m3d.grid.nv
            nvloc = m3d.load.nvloc

            # explicit
            self.nexp = 0
            self.Aexp = []
            if self.use_explicit:
                util.debug(m3d, self, "explicit ...", level=1)
                
                # Name, Count, File
                self.nexp = nexp = self.explicit[1]
                filepath = self.path + self.explicit[2]
                
                # we assume petsc mat aij format in a hdf5 file
                f = util.get_hdf5_file(m3d, self, filepath)
                nrow = f["nrow"].size
                ncol = f["ncol"].size
                nnz = f["nnz"].size
                nnzrow = f["nnzrow"]
                nnzrow_max = numpy.max(nnzrow)
                
                # however, we use index pointers for the petsc call
                indptr = numpy.append(0, numpy.cumsum(nnzrow)).astype("i4")
                colidx = f["colidx"]
                if not nv == nrow == ncol:
                    util.error(m3d, "Dimensions of explicit matrix: {} do not match vector length: {}".format((nrow, ncol), nv))
                    sys.exit(1)
                
                for i in range(nexp):
                    aij = f["aij"][i,...]
                    A = PETSc.Mat()
                    A.create()
                    A.setType("aij")                            # PETSc.Mat.Type.AIJ,
                    A.setSizes(((nvloc, nv), (nvloc, nv)))
                    A.setUp()

                    start, end = A.owner_range
                    A.setPreallocationNNZ((nnzrow_max, nnzrow_max))
                    A.setValuesCSR(
                                   indptr[start:end+1]-indptr[start],
                                   colidx[indptr[start]:indptr[end]],
                                   aij[indptr[start]:indptr[end]])
                    A.assemble()
                    
                    self.Aexp.append(A)

            # implcit
            self.nimp = 0
            self.Aimp = []
            if self.use_implicit:
                util.debug(m3d, self, "implicit ...", level=1)
                
                self.nimp = nimp = self.implicit[1]
                filepath = self.path + self.implicit[2]

                # we assume petsc mat aij format in a hdf5 file
                f = util.get_hdf5_file(m3d, self, filepath)
                nrow = f["nrow"].size
                ncol = f["ncol"].size
                nnz = f["nnz"].size
                nnzrow = f["nnzrow"]
                nnzrow_max = numpy.max(nnzrow)

                # however, we use index pointers for the petsc call
                indptr = numpy.append(0, numpy.cumsum(nnzrow)).astype("i4")
                colidx = f["colidx"]
                if not nv == nrow == ncol:
                    util.error(m3d, "Dimensions of implicit matrix: {} do not match vector length: {}".format((nrow, ncol), nv))
                    sys.exit(1)

                for i in range(nimp):
                    aij = f["aij"][i,...]
                    A = PETSc.Mat()
                    A.create()
                    A.setType("aij")                            # PETSc.Mat.Type.AIJ,
                    A.setSizes(((nvloc, nv), (nvloc, nv)))
                    A.setUp()
                            
                    start, end = A.owner_range
                    A.setPreallocationNNZ((nnzrow_max, nnzrow_max))
                    A.setValuesCSR(
                                   indptr[start:end+1]-indptr[start],
                                   colidx[indptr[start]:indptr[end]],
                                   aij[indptr[start]:indptr[end]])
                    A.assemble()
                    
                    self.Aimp.append(A)


# ----------------------------------------------------------------------------------------
    def __str__(self):

        if not self.use_tmm:
            return "none"

        text = ""
        text = text + "path: {}\n".format(self.path)

        text = text + "explicit:\n"
        if self.use_explicit:
            # Name, Count, File
            text = text + "  {:16} {:5} {:16}\n".format("name", "count", "file")
            explicit = self.explicit

            name = str(explicit[0])
            count = int(explicit[1])
            file = str(explicit[2])
            text = text + "  {:16.16} {:>5d} {:16.16}\n".format(name, count, file)

        else:
            text = text + "  none\n"

        text = text + "implicit:\n"
        if self.use_implicit:
            # Name, Count, File
            text = text + "  {:16} {:5} {:16}\n".format("name", "count", "file")
            implicit = self.implicit

            name = str(implicit[0])
            count = int(implicit[1])
            file = str(implicit[2])
            text = text + "  {:16.16} {:>5d} {:16.16}\n".format(name, count, file)

        else:
            text = text + "  none\n"

        text = text.rstrip()

        return text

## ----------------------------------------------------------------------------------------
#    def apply(self, m3d, Aj, yj):
#
#            m3d.tmm.apply(m3d, Aexpj, yj)





#            nv = m3d.grid.nv.astype("i4")
#            nvloc = m3d.load.nvloc.astype("i4")
#            nvprev = m3d.load.nvprev
#            nv = m3d.grid.nv.astype(PETSc.IntType)
#            nvloc = m3d.load.nvloc.astype(PETSc.IntType)
#                nrow = f["nrow"].size.astype("i4")
#                ncol = f["ncol"].size.astype("i4")
#                nnz = f["nnz"].size.astype("i4")
#                nrow = f["nrow"].size.astype(PETSc.IntType)
#                ncol = f["ncol"].size.astype(PETSc.IntType)
#                nnz = f["nnz"].size.astype(PETSc.IntType)

#                print(nrow.dtype)
#                print(ncol.dtype)
#                print(nnz.dtype)

#                nnzrow = f["nnzrow"][...].astype(PETSc.IntType)
#                nnzrow = f["nnzrow"]
#                rows = numpy.repeat(numpy.arange(nrow), nnzrow)
#                print(rows)
#                print(rows.shape)
#                print(nnzrow.dtype)
#                indptr = numpy.append(0, numpy.cumsum(nnzrow)).astype(PETSc.IntType)
#                colidx = f["colidx"][...].astype(PETSc.IntType)
#                indptr = numpy.append(0, numpy.cumsum(nnzrow))
#                colidx = f["colidx"][...].astype("i4")
#                print(colidx.dtype)

#                    aij = f["aij"][i,...].astype(PETSc.ScalarType)
#        #for ie in nexp:

#        self.nexp = nexp = exp_list[1]
#        self.Aexp = exp_list[2]
#        # !!! read in !!!
#        # !!! read in !!!
#        for i in range(nexp):
#            pass
#        # Aexpj
#
#        self.nimp = nimp = imp_list[1]
#        self.Aimp = imp_list[2]
#        # !!! read in !!!
#        # !!! read in !!!
#        for i in range(nimp):
#            pass
#        # Aimpj

#                    A.setPreallocationNNZ(nnzrow)

##                    A.setSizes(((2, 4), (2, 4)))
#                    A.setSizes(4)
#                    A.setPreallocationNNZ(4)

#                    for j in range(start, end):
###                    for j in range(start, start+100):
##                        print(j)
##                        print(indptr[j])
##                        print(indptr[j+1])
##                        print(colidx[indptr[j]:indptr[j+1]])
#                        A.setValues(j, colidx[indptr[j]:indptr[j+1]], aij[indptr[j]:indptr[j+1]])
##                        A.setValuesLocal(j, colidx[indptr[j]:indptr[j+1]], aij[indptr[j]:indptr[j+1]])

#                    A.setValue(0, 0, 1.)
#                    A.setValues([0,1], [0,1], [1.,2.,3.,4])
#                    A.setValuesCSR([0,2,4,4,4], [0,1,0,1], [1.,2.,3.,4.])
#                    A.setValuesCSR(indptr, colidx, aij)

#                    print(len(indptr[start:end+1]))
#                    print(len(colidx[indptr[start]:indptr[end]]))
#                    print(len(aij[indptr[start]:indptr[end]]))
#                    print(len(indptr[start:end+1]))
#                    sys.stdout.flush()

#                    viewer = PETSc.Viewer().createBinary('matrix-A.dat', 'w')
#                    viewer(A)
#                    A.setValuesIJV(indptr, colidx, aij)

#                    A.view()
#                    print(A.sizes)
#                    sys.exit(1)
