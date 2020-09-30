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
# this is a modified file from
# Davide Scaramuzzas Toolbox OcamCalib
# original filename: findcenter.m

from .calib_data import CalibData
from .calibrate import calibrate
from .libsmop import *


@function
def findcenterUrban(calib_data: CalibData):
    if logical_or(isempty(calib_data.ima_proc), isempty(calib_data.Xp_abs)):
        print('\nNo corner data available. Extract grid corners before calibrating.\n')
        return

    print('\nComputing center coordinates.\n')
    if isempty(calib_data.taylor_order):
        calib_data.taylor_order = calib_data.taylor_order_default

    ## ================ 
    #  added code
    options = optimset('Display', 'off', 'MaxIter', 10000, 'LargeScale', 'off')
    x0 = concat([calib_data.ocam_model.xc, calib_data.ocam_model.yc])
    x0, __, __, __ = fminsearch(errCenterUrban, x0, options, calib_data, nargout=4)
    calib_data.ocam_model.xc = x0[1]
    calib_data.ocam_model.yc = x0[2]
    calib_data.ocam_model.c = 1
    calib_data.ocam_model.d = 0
    calib_data.ocam_model.e = 0
    ## do calibration again
    calib_data.RRfin, calib_data.ocam_model.ss = calibrate(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs,
                                                           calib_data.Yp_abs, calib_data.ocam_model.xc,
                                                           calib_data.ocam_model.yc, calib_data.taylor_order,
                                                           calib_data.ima_proc, nargout=2)
    calib_data.calibrated = 1

    # reproject
    M = concat([calib_data.Xt, calib_data.Yt, zeros(size(calib_data.Xt))])
    allerr, rms = reprojectpoints_advUrban(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc,
                                           calib_data.Xp_abs, calib_data.Yp_abs, M, nargout=2)
    calib_data.rmsAfterCenter = rms
    ss = calib_data.ocam_model.ss
    ss
    figure[3]
    set(3, 'Name', 'Calibration results', 'NumberTitle', 'off')
    subplot(2, 1, 1)
    plot(arange(0, floor(calib_data.width / 2)),
         polyval(concat([ss(arange(end(), 1, - 1))]), concat([arange(0, floor(calib_data.width / 2))])))
    grid('on')
    axis('equal')
    xlabel('Distance \'rho\' from the image center in pixels')
    ylabel('f(rho)')
    title('Forward projection function')

    subplot(2, 1, 2)
    plot(arange(0, floor(calib_data.width / 2)), dot(180 / np.pi, np.arctan2(
        arange(0, floor(calib_data.width / 2)), - polyval(concat([ss(arange(end(), 1, - 1))]), concat(
            [arange(0, floor(calib_data.width / 2))])))) - 90)
    grid('on')
    xlabel('Distance \'rho\' from the image center in pixels')
    ylabel('Degrees')
    title('Angle of optical ray as a function of distance from circle center (pixels)')
    calib_data.calibrated = 1
