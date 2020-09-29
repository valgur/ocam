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


@function
def prova3(x=None, xc=None, yc=None, ss=None, RRfin=None, ima_proc=None, Xp_abs=None, Yp_abs=None, M=None, width=None,
           height=None):
    a = x[1]
    b = x[2]
    c = x[3]
    d = x[4]
    e = x[5]
    ssc = x(arange(6, end()))
    M[:, 3] = 1
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    for i in ima_proc.flat:
        Mc = concat([Mc, dot(RRfin[:, :, i], M.T)])
        Xpp = concat([[Xpp], [Xp_abs[:, :, i]]])
        Ypp = concat([[Ypp], [Yp_abs[:, :, i]]])

    xp1, yp1 = omni3d2pixel(multiply(ss, ssc.T), Mc, width, height, nargout=2)
    xp = dot(xp1, c) + dot(yp1, d) + dot(xc, a)
    yp = dot(xp1, e) + yp1 + dot(yc, b)
    err = sqrt((Xpp - xp.T)**2 + (Ypp - yp.T)**2)
    return err
