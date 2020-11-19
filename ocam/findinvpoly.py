# FINDINVPOLY finds the inverse polynomial specified in the argument.
#   [POL, ERR, N] = FINDINVPOLY(SS, RADIUS, N) finds an approximation of the inverse polynomial specified in OCAM_MODEL.SS.
#   The returned polynomial POL is used in WORLD2CAM_FAST to compute the reprojected point very efficiently.
#   
#   SS is the polynomial which describe the mirrror/lens model.
#   RADIUS is the radius (pixels) of the omnidirectional picture.
#   ERR is the error (pixel) that you commit in using the returned
#   polynomial instead of the inverse SS. N is searched so that
#   that ERR is < 0.01 pixels.
#
#   Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
#   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org

import numpy as np
import cv2


def findinvpoly(ss, radius, N=None, tol=0.01):
    if N is None:
        maxerr = np.inf
        N = 1
        while maxerr > tol:
            N += 1
            pol, err, _ = findinvpoly2(ss, radius, N)
            maxerr = np.max(err)
    else:
        pol, err, N = findinvpoly2(ss, radius, N)

    return pol, err, N


def findinvpoly2(ss, radius, N):
    theta = np.arange(-np.pi / 2, 1.2, 0.01)
    r = invFUN(ss, theta)
    valid = r < radius
    theta = theta[valid]
    r = r[valid]
    pol = np.polyfit(theta, r, N)
    err = np.abs(r - np.polyval(pol, theta))
    return pol, err, N


def invFUN(ss, theta):
    m = np.tan(theta)
    r = np.zeros(len(m))
    poly_coef = ss[::-1]
    poly_coef_tmp = poly_coef.copy()
    for i in range(len(m)):
        poly_coef_tmp[-2] = poly_coef[-2] - m[i]
        rho_tmp = cv2.solvePoly(poly_coef_tmp[::-1], maxIters=50)[1][:, 0]
        res = rho_tmp[(rho_tmp[:, 0] > 0) & (np.abs(rho_tmp[:, 1]) < 1e-10), 0]
        r[i] = np.min(res) if len(res) > 0 else np.inf
    return r
