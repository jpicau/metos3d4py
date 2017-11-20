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

class TimeStep:
    """
        TimeStep class
        
        TimeStep context
        
        Attributes:
            t0          start
            nt          count
            dt          step

        """

# ----------------------------------------------------------------------------------------

    def init(self, comm, conf):
        
        self.t0 = conf.dict["Time step"]["Start"]
        self.nt = conf.dict["Time step"]["Count"]
        self.dt = conf.dict["Time step"]["Step"]

# ----------------------------------------------------------------------------------------

    def __str__(self):
        return "TimeStep:\n  t0: {}\n  nt: {}\n  dt: {}".format(self.t0, self.nt, self.dt)
        


