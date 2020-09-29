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

from .libsmop import *

from .omni3d2pixel import omni3d2pixel
from .rodrigues import rodrigues


@function
def prova(x=None, ss=None, int_par=None, Xp_abs=None, Yp_abs=None, M=None, width=None, height=None):
    c = int_par[1]
    d = int_par[2]
    e = int_par[3]
    xc = int_par[4]
    yc = int_par[5]
    num_points = size(M, 1)
    R = rodrigues(concat([x[1], x[2], x[3]]))
    T = concat([x[4], x[5], x[6]]).T
    Mc = dot(R, M.T) + dot(T, ones(1, num_points))
    xp1, yp1 = omni3d2pixel(ss, Mc, width, height, nargout=2)
    xp = dot(xp1, c) + dot(yp1, d) + xc
    yp = dot(xp1, e) + yp1 + yc
    err = sqrt((Xp_abs - xp.T)**2 + (Yp_abs - yp.T)**2)
    return err
