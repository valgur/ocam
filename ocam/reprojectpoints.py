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
def reprojectpoints(calib_data: CalibData):
    m = copy([])
    xx = copy([])
    err = copy([])
    stderr = copy([])
    rhos = copy([])
    num_points = size(calib_data.Xp_abs, 1)
    MSE = 0
    counterr = 0
    for i in calib_data.ima_proc.flat:
        counterr += 1
        xx = dot(calib_data.RRfin[:, :, i],
                 concat([[calib_data.Xt.T], [calib_data.Yt.T], [ones(size(calib_data.Xt.T))]]))
        Xp_reprojected, Yp_reprojected = omni3d2pixel(calib_data.ocam_model.ss, xx, calib_data.width,
                                                      calib_data.height, nargout=2)
        stt = sqrt((calib_data.Xp_abs[:, :, i] - calib_data.ocam_model.xc - Xp_reprojected.T)**2 + (
                calib_data.Yp_abs[:, :, i] - calib_data.ocam_model.yc - Yp_reprojected.T)**2)
        err[counterr] = mean(stt)
        stderr[counterr] = std(stt)
        MSE += sum((calib_data.Xp_abs[:, :, i] - calib_data.ocam_model.xc - Xp_reprojected.T)**2 + (
                calib_data.Yp_abs[:, :, i] - calib_data.ocam_model.yc - Yp_reprojected.T)**2)

    print('\n Average reprojection error computed for each chessboard [pixels]:\n')
    for i in arange(1, length(err)).flat:
        fprintf(' %3.2f %c %3.2f\n', err[i], 177, stderr[i])

    # err'
    print('\n Average error [pixels]\n\n %f\n', mean(err), end='')
    print('\n Sum of squared errors\n\n %f\n', MSE, end='')
