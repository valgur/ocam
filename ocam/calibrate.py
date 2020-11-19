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
import numpy.linalg as la
from scipy.special import factorial

from ocam.lsqlin import lsqlin


def calibrate_linear(Xt, Yt, Xp_abs, Yp_abs, xc, yc, taylor_order, ima_proc):
    Yp = Yp_abs - yc
    Xp = Xp_abs - xc
    n_img, n_corners = Xp_abs.shape
    RRfin = np.zeros((n_img, 3, 3))
    for kk in ima_proc:
        Ypt = Yp[kk, :]
        Xpt = Xp[kk, :]
        # Building the matrix
        A = np.c_[Xt * Ypt, Yt * Ypt, -Xt * Xpt, -Yt * Xpt, Ypt, -Xpt]
        U, S, Vt = la.svd(A, full_matrices=False)
        # Solving for computing the scale factor Lambda
        R11, R12, R21, R22, T1, T2 = Vt[-1, :]

        AA = (R11 * R12 + R21 * R22)**2
        BB = R11**2 + R21**2
        CC = R12**2 + R22**2
        R32_2 = np.roots([1, CC - BB, -AA])
        R32_2 = R32_2[R32_2 >= 0]

        if R32_2 != 0:  # <= 0.00000001 * (R12^2+R22^2) %that is like R32==0
            R32 = np.sqrt(R32_2)
            R32 = np.c_[R32, -R32].flatten()
            R31 = -(R11 * R12 + R21 * R22) / R32
        else:
            R32 = np.sqrt(R32_2)
            R32 = np.c_[R32, R32, -R32, -R32].flatten()
            R31 = np.repeat(np.sqrt(CC - BB), len(R32_2))
            R31 = np.c_[R31, -R31, R31, -R31].flatten()

        RR = np.zeros((len(R32), 3, 3))
        for i in range(len(R32)):
            Lb = 1. / np.sqrt(R11**2 + R21**2 + R31[i]**2)
            RR[i] = Lb * np.array([
                [R11, R12, T1],
                [R21, R22, T2],
                [R31[i], R32[i], 0]
            ])
        RR = np.r_[RR, -RR]

        T = RR[:, :2, 2]
        minRR_ind = np.argmin(la.norm(T - np.array([Xpt[0], Ypt[0]]), axis=1))
        RR1 = RR[np.all(np.sign(T) == np.sign(T[minRR_ind]), axis=1)]
        if len(RR1) == 0:
            return None, None
        nm = find_valid_RR(RR1, Xt, Yt, Xpt, Ypt)
        RRfin[kk] = RR1[nm]

    # find mirror and translation parameters
    RRfin, ss = omni_find_parameters_fun(Xt, Yt, Xp_abs, Yp_abs, xc, yc, RRfin, taylor_order, ima_proc)
    return RRfin, ss


def find_valid_RR(RR, Xt, Yt, Xpt, Ypt):
    rho2 = Xpt**2 + Ypt**2
    rho = np.sqrt(rho2)
    for i in range(len(RR)):
        (R11, R12, T1,
         R21, R22, T2,
         R31, R32, _) = RR[i].flat
        MA = R21 * Xt + R22 * Yt + T2
        MB = Ypt * (R31 * Xt + R32 * Yt)
        MC = R11 * Xt + R12 * Yt + T1
        MD = Xpt * (R31 * Xt + R32 * Yt)
        PP1 = np.r_[
            np.c_[MA, MA * rho, MA * rho2],
            np.c_[MC, MC * rho, MC * rho2]
        ]
        PP = np.c_[PP1, np.r_[-Ypt, -Xpt]]
        QQ = np.r_[MB, MD]
        s = np.squeeze(la.pinv(PP) @ QQ)
        ss = s[:3]
        if ss[2] >= 0:
            return i
    return None


def omni_find_parameters_fun(Xt, Yt, Xp_abs, Yp_abs, xc, yc, RRfin, taylor_order, ima_proc,
                             ensure_convexity=True):
    Yp = Yp_abs - yc
    Xp = Xp_abs - xc
    n_img, n_corners = Xp_abs.shape
    PP = np.zeros((n_img, 2 * n_corners, taylor_order + n_img))
    QQ = np.zeros((n_img, 2, n_corners))
    for i, kk in enumerate(ima_proc):
        (R11, R12, T1,
         R21, R22, T2,
         R31, R32, _) = RRfin[kk].flat
        Xpt = Xp[kk]
        Ypt = Yp[kk]
        MA = R21 * Xt + R22 * Yt + T2
        MB = Ypt * (R31 * Xt + R32 * Yt)
        MC = R11 * Xt + R12 * Yt + T1
        MD = Xpt * (R31 * Xt + R32 * Yt)
        rho = np.zeros((n_corners, taylor_order))
        rho[:, 1:] = np.sqrt(Xpt**2 + Ypt**2)[None].T**np.arange(2, taylor_order + 1)
        PP[i, :, :taylor_order] = np.r_[MA[None].T * rho, MC[None].T * rho]
        PP[i, :, 0] = np.r_[MA, MC]
        PP[i, :, taylor_order + i] = np.r_[-Ypt, -Xpt]
        QQ[i, 0] = MB
        QQ[i, 1] = MD
    PP = PP.reshape(-1, PP.shape[-1])
    QQ = QQ.reshape(-1, 1)

    if ensure_convexity:
        # Ensure that the resulting polynomial is convex by adding a non-negative constraint on its second derivative.
        # This guarantees that the inverse of the polynomial can be mapped nicely one-to-one to another polynomial.
        rho = np.arange(int(np.ceil(Xp_abs.max())))
        mono, mono_const = construct_convexity_constraint(n_img, rho, taylor_order)
        solution = lsqlin(PP, QQ, A=mono, b=mono_const)
        assert solution['status'] == 'optimal'
        s = np.squeeze(solution['x'])
    else:
        s = np.squeeze(la.pinv(PP) @ QQ)

    ss = s[:taylor_order]
    for i, kk in enumerate(ima_proc):
        RRfin[kk, 2, 2] = s[taylor_order + i]
    ss = np.r_[ss[0], 0, ss[1:]]
    return RRfin, ss


def construct_convexity_constraint(n_img, rho, taylor_order):
    mono = np.empty((0, taylor_order))
    diff_order = 2
    for diff_k in range(1, min(diff_order, taylor_order) + 1):
        mono1 = []
        for i in range(diff_k):
            mono1.append(np.zeros_like(rho))
        mono1.append(np.full_like(rho, factorial(diff_k)))
        for i in range(diff_k + 1, taylor_order + 1):
            mono1.append((factorial(i) / factorial(i - diff_k)) * rho**(i - diff_k))
        del mono1[1]
        mono = np.r_[mono, np.array(mono1).T]
    if len(mono) == 0:
        return None, None
    mono = np.c_[-mono, np.zeros((len(mono), n_img))]
    mono_const = np.zeros((len(mono), 1))
    return mono, mono_const
