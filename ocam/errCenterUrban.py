#     Steffen Urban email: steffen.urban@kit.edu
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
#
# 04.03.2014 by Steffen Urban
# error function for the center of distortion search

from math import isinf

from .calib_data import CalibData
from .calibrate import calibrate
from .libsmop import *
from .omni3d2pixel import omni3d2pixel


@function
def errCenterUrban(x, calib_data: CalibData):
    error = 0
    xc = x[1]
    yc = x[2]
    # call calibration function
    RRfin, ss = calibrate(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs, calib_data.Yp_abs, xc, yc,
                          calib_data.taylor_order, calib_data.ima_proc, nargout=2)
    lauf = 1
    M = concat([calib_data.Xt, calib_data.Yt, ones(size(calib_data.Xt))])
    for i in arange(1, size(RRfin, 3)).flat:
        # if calibration was not possible add a high value
        # to penalize the minimization away from that point
        if calib_data.RRfin[:, :, i] == 0:
            error += sum(dot(ones(length(calib_data.Xp_abs), 1),
                             sqrt((calib_data.width / 2)**2 + (calib_data.height / 2)**2)))
        else:
            Mc = dot(RRfin[:, :, i], M.T)
            Xpp = calib_data.Xp_abs[:, :, i]
            Ypp = calib_data.Yp_abs[:, :, i]
            xp1, yp1 = omni3d2pixel(ss, Mc, calib_data.width, calib_data.height, nargout=2)
            if (logical_or(isinf(xp1), isinf(yp1))):
                error += sum(dot(ones(length(calib_data.Xp_abs), 1),
                                 sqrt((calib_data.width / 2)**2 + (calib_data.height / 2)**2)))
            else:
                xp = xp1 + xc
                yp = yp1 + yc
                lauf += length(Xpp)
                error += sum((Xpp - xp.T)**2) + sum((Ypp - yp.T)**2)

    error = sqrt(error / lauf)
    return error
