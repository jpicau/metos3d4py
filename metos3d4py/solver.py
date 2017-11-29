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

        self.yl = []
        for y in m3d.tracer.y0:
            self.yl.append(y.duplicate())

        # debug
        util.debug(m3d, self, self, level=1)


### ----------------------------------------------------------------------------------------
##    def solve_empty(self, m3d):
##        util.debug(m3d, self, "Testing empty solver. Just returns.")
##        return

# ----------------------------------------------------------------------------------------
    def solve(self, m3d):
        util.debug(m3d, self, "Solve ...", level=1)

        tracer = m3d.tracer
        time = m3d.time

        ny = tracer.ny
        y0 = tracer.y0
        
        nl = self.nl
        yl = self.yl
        
        i = 0
#        print("rank: {}, i: {}, y[i]: {}, y0[i]: {}".format(m3d.rank, i, yl[i], y0[i]))
        util.copy_vector_list(y0, yl)
#        print("rank: {}, i: {}, y[i]: {}, y0[i]: {}".format(m3d.rank, i, yl[i], y0[i]))
#        for i in range(ny):
#            print("rank: {}, i: {}, y[i]: {}, y0[i]: {}".format(m3d.rank, i, yl[i], y0[i]))
##            yl[i] = y0[i]
#            y0[i].copy(yl[i])
##            yl[i].copy(y0[i])
#            print("rank: {}, i: {}, y[i]: {}, y0[i]: {}".format(m3d.rank, i, yl[i], y0[i]))
#        import sys
#        sys.exit(1)

        for i in range(nl):
            util.debug(m3d, self, "{}".format(i), level=1)
            time.step(m3d, yl)

# ----------------------------------------------------------------------------------------
    def __str__(self):
        text = ""
        text = text + "   type: {}\n".format(self.type)
        text = text + "  count: {}\n".format(self.nl)
        text = text + "monitor: {}".format(self.mon)
        return text






