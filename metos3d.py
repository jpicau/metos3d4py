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

"""
    Metos3D
    =======
    
    Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
    
    [Piwonski and Slawig, 2016]
    https://www.geosci-model-dev.net/9/3729/2016/
    
"""

import sys
import metos3d4py as m3d

if __name__ == "__main__":
    
    sim = m3d.Metos3D(sys.argv)
    sim.run()


