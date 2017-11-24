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
        return "Time:\n  t0: {}\n  nt: {}\n  dt: {}".format(self.t0, self.nt, self.dt)

## ----------------------------------------------------------------------------------------
#    def get(self):
#        return self.t0, self.nt, self.dt

# ----------------------------------------------------------------------------------------
    def init(self, m3d):
        
#        self.bj = bj
#        self.dj = dj
#        self.yj = yj
#        self.qj = qj

        comm = m3d.comm
        conf = m3d.config["Time"]

        self.t0 = conf["Start"]
        self.nt = conf["Count"]
        self.dt = conf["Step"]

        m3d.time = self




