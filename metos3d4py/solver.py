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

from metos3d4py import util

class Solver:
    """
        Solver class
        
        Solver context
        
        Attributes:
            nl          model years
            yl          ...

        """

# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):

        config = util.get_key(m3d, self, m3d.config, "Solver", dict)

        self.type = util.get_key(m3d, self, config, "Type", str)
        self.nl = util.get_key(m3d, self, config, "Count", int)
        self.mon = util.get_key(m3d, self, config, "Monitor", bool)

#        self.yl = []

        # debug
        util.debug(m3d, self, self, level=1)


### ----------------------------------------------------------------------------------------
##    def solve_empty(self, m3d):
##        util.debug(m3d, self, "Testing empty solver. Just returns.")
##        return
#
## ----------------------------------------------------------------------------------------
#    def solve(self, m3d):
#        comm = m3d.comm
#        time = m3d.time
#        util.debug(m3d, self, "Solve ...")
#
#        nl = self.nl
#        yl = self.yl
##        for i in range(nl):

# ----------------------------------------------------------------------------------------
    def __str__(self):
        text = ""
        text = text + "   type: {}\n".format(self.type)
        text = text + "  count: {}\n".format(self.nl)
        text = text + "monitor: {}".format(self.mon)
        return text






