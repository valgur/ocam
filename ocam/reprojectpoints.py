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

import numpy as np

from .calib_data import CalibData
from .world2cam import omni3d2pixel


def reprojectpoints(calib_data: CalibData):
    M = np.c_[calib_data.Xt, calib_data.Yt, np.ones(len(calib_data.Xt))]
    return reprojectpoints_adv(
        calib_data.ocam_model,
        calib_data.RRfin,
        calib_data.ima_proc,
        calib_data.Xp_abs,
        calib_data.Yp_abs,
        M
    )


def reprojectpoints_adv(ocam_model, RRfin, ima_proc, Xp_abs, Yp_abs, M):
    n_img = len(ima_proc)
    err = np.zeros(n_img)
    stderr = np.zeros(n_img)
    mse = 0
    ss = ocam_model.ss
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    xc = ocam_model.xc
    yc = ocam_model.yc
    xx = RRfin @ M.T
    for i, kk in enumerate(ima_proc):
        Xp_reprojected, Yp_reprojected = omni3d2pixel(ss, xx[kk])
        Xp_reprojected = Xp_reprojected * c + Yp_reprojected * d + xc
        Yp_reprojected = Xp_reprojected * e + Yp_reprojected + yc
        sqerr = (Xp_abs[kk] - Xp_reprojected)**2 + (Yp_abs[kk] - Yp_reprojected)**2
        stt = np.sqrt(sqerr)
        err[i] = stt.mean()
        stderr[i] = stt.std()
        mse += sqerr.sum()

    print('\n Average reprojection error computed for each chessboard [pixels]:\n')
    for kk in range(len(err)):
        print(f' {err[kk]:3.2f} {177:c} {stderr[kk]:3.2f}')

    print('\n Average error [pixels]\n ', err.mean())
    print('\n Sum of squared errors\n ', mse)

    return err, stderr, mse
