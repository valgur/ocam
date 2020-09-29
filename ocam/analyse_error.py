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

from .calib_data import CalibData
from .libsmop import *
from .omni3d2pixel import omni3d2pixel


@function
def analyse_error(calib_data: CalibData):
    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    figure[5]
    set(5, 'Name', 'Analyse error', 'NumberTitle', 'off')
    zoom('on')
    colors = 'brgkcm'
    m = copy([])
    xx = copy([])
    err = copy([])
    stderr = copy([])
    rhos = copy([])
    num_points = size(calib_data.Xp_abs, 1)
    MSE = 0
    count = 0
    if logical_and(logical_and(isempty(calib_data.ocam_model.c), isempty(calib_data.ocam_model.d)),
                   isempty(calib_data.ocam_model.e)):
        calib_data.ocam_model.c = 1
        calib_data.ocam_model.d = 0
        calib_data.ocam_model.e = 0

    for i in calib_data.ima_proc.flat:
        count += 1
        xx = dot(calib_data.RRfin[:, :, i],
                 concat([[calib_data.Xt.T], [calib_data.Yt.T], [ones(size(calib_data.Xt.T))]]))
        xp1, yp1 = omni3d2pixel(calib_data.ocam_model.ss, xx, calib_data.ocam_model.width, calib_data.ocam_model.height,
                                nargout=2)
        xp = dot(xp1, calib_data.ocam_model.c) + dot(yp1, calib_data.ocam_model.d) + calib_data.ocam_model.xc
        yp = dot(xp1, calib_data.ocam_model.e) + yp1 + calib_data.ocam_model.yc
        sqerr = (calib_data.Xp_abs[:, :, i] - xp.T)**2 + (
                calib_data.Yp_abs[:, :, i] - yp.T)**2
        err[count] = np.mean(sqrt(sqerr))
        stderr[count] = np.std(sqrt(sqerr))
        MSE += sum(sqerr)
        plot(calib_data.Xp_abs[:, :, i] - xp.T, calib_data.Yp_abs[:, :, i] - yp.T,
             concat([str(colors(rem(i - 1, 6) + 1)), '+']))
        hold('on')

    hold('off')
    grid('on')
    print('\n Average reprojection error computed for each chessboard [pixels]:\n')
    for i in arange(1, length(err)).flat:
        fprintf(' %3.2f %c %3.2f\n', err(i), 177, stderr(i))

    print('\n Average error [pixels]\n\n %f\n', np.mean(err), end='')
    print('\n Sum of squared errors\n\n %f\n', MSE, end='')
    ss = calib_data.ocam_model.ss
    ss
