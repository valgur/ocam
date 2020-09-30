###########################################################################  
#   Copyright (C) 2006 DAVIDE SCARAMUZZA
#   
#   Author: Davide Scaramuzza - email: davsca@tiscali.it
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#   USA
############################################################################

from .FUNrho import FUNrho
from .libsmop import *


# Find the other parameters
@function
def omni_find_extrs_parameters(ss=None, xc=None, yc=None, ima_proc=None, Xp_abs=None, Yp_abs=None, Xt=None, Yt=None,
                               RRfin=None):
    varargin = omni_find_extrs_parameters.varargin
    nargin = omni_find_extrs_parameters.nargin

    Xp = Xp_abs - xc
    Yp = Yp_abs - yc
    M1 = copy([])
    M2 = copy([])
    M3 = copy([])
    Rt = copy([])
    MM = copy([])
    RRdef = copy([])
    for i in ima_proc.flat:
        Xpt = Xp[:, :, i]
        Ypt = Yp[:, :, i]
        rhot = sqrt(Xpt**2 + Ypt**2)
        M1 = concat(
            [zeros(size(Xt)), zeros(size(Xt)), multiply(- FUNrho(ss, rhot), Xt), multiply(- FUNrho(ss, rhot), Yt),
             multiply(Ypt, Xt), multiply(Ypt, Yt), zeros(size(Xt)), - FUNrho(ss, rhot), Ypt])
        M2 = concat([multiply(FUNrho(ss, rhot), Xt), multiply(FUNrho(ss, rhot), Yt), zeros(size(Xt)), zeros(size(Xt)),
                     multiply(- Xpt, Xt), multiply(- Xpt, Yt), FUNrho(ss, rhot), zeros(size(Xt)), - Xpt])
        M3 = concat([multiply(- Ypt, Xt), multiply(- Ypt, Yt), multiply(Xpt, Xt), multiply(Xpt, Yt), zeros(size(Xt)),
                     zeros(size(Xt)), - Ypt, Xpt, zeros(size(Xt))])
        MM = concat([[M1], [M2], [M3]])
        U, S, V = svd(MM)
        res = V[:, end()]
        Rt = reshape(res(arange(1, 6)), 2, 3).T
        scalefact = sqrt(abs(dot(norm(Rt[:, 1]), norm(Rt[:, 2]))))
        #    keyboard;
        Rt = concat([Rt, cross(Rt[:, 1], Rt[:, 2])])
        U2, S2, V2 = svd(Rt)
        Rt = dot(U2, V2.T)
        Rt[:, 3] = res(arange(7, end())) / scalefact
        Rt = dot(Rt, sign(dot(Rt[1, 3], RRfin(1, 3, i))))
        RRdef[:, :, i] = Rt
        # pause

    RRtmp = copy(RRfin)
    RRfin = copy(RRdef)

    return RRfin
