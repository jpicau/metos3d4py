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

from metos3d4py.util import util

class Solver:
    """
        Solver class
        
        Solver context
        
        Attributes:
            nl          model years
            yl          ...

        """

# ----------------------------------------------------------------------------------------
    def solve(self, m3d):
        comm = m3d.comm
        time = m3d.time
        util._print(comm, "Solver: Solve ...")
    
        nl = self.nl
        yl = self.yl

# ----------------------------------------------------------------------------------------
    def __str__(self):
        return "Solver:\n  nl: {}".format(self.nl)

# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        import sys
        # check solver key
        comm = m3d.comm
        try:
            conf = m3d.conf["Solver"]
        except Exception as e:
            util._print_error(comm, "Solver init: Cannot retrieve '{}' key from configuration.".format("Solver"))
            util._print_usage(comm)
            sys.exit(1)
        
        self.nl = conf["Count"]
        self.yl = []

        m3d.solver = self




