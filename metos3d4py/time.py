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

class Time:
    """
        Time step context
        
        Attributes:
            t0          start
            nt          count
            dt          step

        """

# ----------------------------------------------------------------------------------------
    def __str__(self):
        text = ""
        text = text + "start: {}\n".format(self.t0)
        text = text + "count: {}\n".format(self.nt)
        text = text + " step: {}".format(self.dt)
        return text

# ----------------------------------------------------------------------------------------
    def __init__(self, m3d):
        
        config = util.get_key(m3d, self, m3d.config, "Time", dict)

        self.t0 = util.get_key(m3d, self, config, "Start", float)
        self.nt = util.get_key(m3d, self, config, "Count", int)
        self.dt = util.get_key(m3d, self, config, "Step", float)

#        self.bj = bj
#        self.dj = dj
#        self.yj = yj
#        self.yexpj = yexpj
#        self.qj = qj

        # debug
        util.debug(m3d, self, self, level=1)








