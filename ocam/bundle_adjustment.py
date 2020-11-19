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
import cv2
import numpy as np
from scipy.optimize import least_squares

from .calib_data import CalibData
from .reprojectpoints import omni3d2pixel


def bundle_adjustment(calib_data, robust, verbose=True, optim_args=None):
    if robust:
        print('Starting robust non-linear refinement')
    else:
        print('Starting non-linear refinement')

    if calib_data.n_imgs == 0 or not calib_data.calibrated:
        print('\nNo linear estimate available. You must first calibrate your camera.\nClick on "Calibration"\n')
        return

    M = np.c_[calib_data.Xt, calib_data.Yt, np.zeros(len(calib_data.Xt))]
    ss0 = calib_data.ocam_model.ss
    x0 = [1, 1, 1, 0, 0] + [1] * len(calib_data.ocam_model.ss)
    offset = 6 + calib_data.taylor_order
    for kk in calib_data.ima_proc:
        R = calib_data.RRfin[kk].copy()
        R[:, 2] = np.cross(R[:, 0], R[:, 1])
        r = cv2.Rodrigues(R)[0].ravel()
        t = calib_data.RRfin[kk, :, 2]
        x0 += [r[0], r[1], r[2], t[0], t[1], t[2]]
    x0 = np.array(x0)

    default_args = dict(xtol=1e-05, ftol=0.0001, max_nfev=1000, verbose=2 if verbose else 0)
    optim_args = {**default_args, **(optim_args or {})}
    ba_result = least_squares(bundle_err, x0, args=[calib_data, M, robust], **optim_args)
    x_opt = ba_result['x']
    RT = x_opt[offset:].reshape(-1, 2, 3)
    RRfinOpt = np.zeros_like(calib_data.RRfin)
    for i, kk in enumerate(calib_data.ima_proc):
        RRfinOpt[kk], _ = cv2.Rodrigues(RT[i, 0])
        RRfinOpt[kk, :, 2] = RT[i, 1]

    ssc = x_opt[5:offset]
    calib_data.ocam_model.ss = ss0 * ssc
    calib_data.ocam_model.xc *= x_opt[0]
    calib_data.ocam_model.yc *= x_opt[1]
    calib_data.ocam_model.c = x_opt[2]
    calib_data.ocam_model.d = x_opt[3]
    calib_data.ocam_model.e = x_opt[4]
    calib_data.RRfin = RRfinOpt
    calib_data.ba_result = ba_result
    calib_data.optimized = True


def bundle_err(x: np.ndarray, calib_data: CalibData, M: np.ndarray, robust: bool):
    a, b, c, d, e = x[:5]
    offset = 6 + calib_data.taylor_order
    ssc = x[5:offset]
    n_img = len(calib_data.ima_proc)
    n_corners = len(M)
    Mc = np.empty((3, n_img, n_corners))
    RT = x[offset:].reshape(-1, 2, 3)
    for i in calib_data.ima_proc:
        R, _ = cv2.Rodrigues(RT[i, 0])
        T = RT[i, 1]
        Mc[:, i, :] = R @ M.T + T[None].T
    Mc = Mc.reshape(3, -1)

    xp1, yp1 = omni3d2pixel(calib_data.ocam_model.ss * ssc, Mc)
    xp = calib_data.ocam_model.xc * a + xp1 * c + yp1 * d
    yp = calib_data.ocam_model.yc * b + xp1 * e + yp1
    errx = calib_data.Xp_abs.ravel() - xp
    erry = calib_data.Yp_abs.ravel() - yp
    errW = np.r_[errx, erry]
    if robust:
        w = huber_weight(errW)
        calib_data.weights = w
        errW = np.sqrt(w) * errW

    errW /= n_img
    return errW


def huber_weight(v, k=1):
    a = np.abs(v) <= k
    b = ~a
    weight = np.empty_like(v)
    weight[a] = v[a] / v[a]
    weight[b] = k / np.abs(v[b])
    weight = weight.T
    return weight
