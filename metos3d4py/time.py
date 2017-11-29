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

        tracer = m3d.tracer
        ny = tracer.ny
        y0 = tracer.y0
        
        bgc = m3d.bgc
        b = bgc.b
        nb = bgc.nb
        nbi = bgc.nbi
        bj = []
        for bi in b:
            bj.append(bi[0].duplicate())
        self.bj = bj
        
        d = bgc.d
        nd = bgc.nd
        ndi = bgc.ndi
        dj = []
        for di in d:
            dj.append(di[0].duplicate())
        self.dj = dj

        yj = []
        qj = []
        for yi in y0:
            yj.append(yi.duplicate())
            qj.append(yi.duplicate())
        self.yj = yj
        self.qj = qj
#        self.yexpj = yexpj

        # debug
        util.debug(m3d, self, self, level=1)

# ----------------------------------------------------------------------------------------
    def step(self, m3d, yl):
        util.debug(m3d, self, "time step ...", level=1)

        tracer = m3d.tracer
        ny = tracer.ny
        y0 = tracer.y0

        bgc = m3d.bgc
        u = bgc.u
        
        tmm = m3d.tmm
        
        t0 = self.t0
        dt = self.dt
        nt = self.nt
        yj = self.yj
        qj = self.qj
        bj = self.bj
        dj = self.dj

        for i in range(ny):
            yj[i] = yl[i]

        for j in range(nt):
            util.debug(m3d, self, "{}".format(i), level=1)
            tj = t0 + j*dt

            if bgc.use_bgc:
                util.debug(m3d, self, "use_bgc", level=1)
                
                if bgc.use_boundary:
                    b = bgc.b
                    nb = bgc.nb
                    nbi = bgc.nbi
                    # interpolate
                    for ib in range(nb):
                        alpha, ialpha, beta, ibeta = util.interpolate(nbi[ib], tj)
                        bj[ib] = alpha * b[ib][ialpha] + beta * b[ib][ibeta]
            
                if bgc.use_domain:
                    d = bgc.d
                    nd = bgc.nd
                    ndi = bgc.ndi
                    # interpolate
                    for id in range(nd):
                        alpha, ialpha, beta, ibeta = util.interpolate(ndi[id], tj)
                        dj[id] = alpha * d[id][ialpha] + beta * d[id][ibeta]

                m3d.bgc.evaluate(m3d, dt, qj, tj, yj, u, bj, dj)

            if tmm.use_tmm:
                util.debug(m3d, self, "use_tmm", level=1)
                
                if tmm.use_explicit:
                    # interpolate
                    
                    # apply explicit
                    m3d.tmm.apply_explicit(m3d, tj, yj, qj)
                
                if tmm.use_implicit:
                    # interpolate
                    
                    # apply implicit
                    m3d.tmm.apply_implicit(m3d, tj, yj, qj)







