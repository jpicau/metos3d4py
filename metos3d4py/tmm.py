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

        # Name, Count, File
        self.nexp = explicit[1]

# ----------------------------------------------------------------------------------------
    def init_implicit(self, m3d):

        self.implicit = implicit = self.config.get("Implicit")
        if implicit is None:
            self.use_implicit = False
            return
        else:
            self.use_implicit = True

        # Name, Count, File
        self.nimp = implicit[1]

# ----------------------------------------------------------------------------------------
    def set(self):
        util.debug(m3d, self, "TMM set ...", level=1)

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



