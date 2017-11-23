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

    def init(self, m3d):
        
        comm = m3d.comm
        conf_tmm = m3d.config["TMM"]
    
        self.path = path = conf_tmm["Path"]
        
        # Name, Count, File
        self.exp_list = exp_list = conf_tmm["Explicit"]
        self.imp_list = imp_list = conf_tmm["Implicit"]
    
        self.nexp = nexp = exp_list[1]
        self.Aexp = exp_list[2]
        # !!! read in !!!
        # !!! read in !!!
        for i in range(nexp):
            pass
        # Aexpj

        self.nimp = nimp = imp_list[1]
        self.Aimp = imp_list[2]
        # !!! read in !!!
        # !!! read in !!!
        for i in range(nimp):
            pass
        # Aimpj

        m3d.tmm = self

# ----------------------------------------------------------------------------------------

    def __str__(self):
        return "TMM:\n  nexp: {}\n  nimp: {}".format(self.nexp, self.nimp)
        


