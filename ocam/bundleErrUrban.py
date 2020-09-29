#     Steffen Urban
#     Copyright (C) 2014  Steffen Urban
# 
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 04.03.2014 by Steffen Urban

from .calib_data import CalibData
from .libsmop import *
from .omni3d2pixel import omni3d2pixel
from .rodrigues import rodrigues


@function
def bundleErrUrban(x, calib_data: CalibData, M, robust):
    global weights
    a = x[1]
    b = x[2]
    c = x[3]
    d = x[4]
    e = x[5]
    offset = calib_data.taylor_order + 6
    ssc = x(arange(6, offset))
    num_points = size(M, 1)
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    lauf = 0
    for i in calib_data.ima_proc.flat:
        R = rodrigues(concat([x(offset + 1 + lauf), x(offset + 2 + lauf), x(offset + 3 + lauf)]))
        T = concat([x(offset + 4 + lauf), x(offset + 5 + lauf), x(offset + 6 + lauf)]).T
        Mc = concat([Mc, dot(R, M.T) + dot(T, ones(1, num_points))])
        Xpp = concat([[Xpp], [calib_data.Xp_abs[:, :, i]]])
        Ypp = concat([[Ypp], [calib_data.Yp_abs[:, :, i]]])
        lauf += 6

    xp1, yp1 = omni3d2pixel(multiply(calib_data.ocam_model.ss, ssc.T), Mc, calib_data.ocam_model.width,
                            calib_data.ocam_model.height, nargout=2)
    xp = dot(xp1, c) + dot(yp1, d) + dot(calib_data.ocam_model.xc, a)
    yp = dot(xp1, e) + yp1 + dot(calib_data.ocam_model.yc, b)
    errx = Xpp - xp.T
    erry = Ypp - yp.T
    if (logical_not(robust)):
        errW = concat([errx, erry])
    else:
        errn = concat([errx, erry])
        w = huberWeight(errn)
        weights = copy(w)
        errW = multiply(sqrt(w), errn)

    return errW


@function
def huberWeight(v):
    k = 1
    a = abs(v) <= k
    b = abs(v) > k
    weight = copy([])
    weight[find(a)] = v(a) / v(a)
    weight[find(b)] = k / abs(v(b))
    weight = weight.T
    return weight
