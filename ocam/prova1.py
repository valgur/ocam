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
def prova1(ss=None, RRfin=None, Xp=None, Yp=None, M=None):
    num_points = size(M, 1)
    Nframes = size(RRfin, 3)
    M[:, 3] = 1
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    for i in arange(1, Nframes).flat:
        Mc = concat([Mc, dot(RRfin[:, :, i], M.T)])
        Xpp = concat([[Xpp], [Xp[:, :, i]]])
        Ypp = concat([[Ypp], [Yp[:, :, i]]])

    xp, yp = omni3d2pixel(ss, Mc, nargout=2)
    err = sqrt((Xpp - xp.T)**2 + (Ypp - yp.T)**2)
    return err
