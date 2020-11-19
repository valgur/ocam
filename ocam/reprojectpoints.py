import cv2
import numpy as np
from numpy import linalg as la

from .calib_data import CalibData, OCamModel

eps = np.finfo(float).eps


def world2cam(M, ocam_model: OCamModel):
    """Projects a 3D point on to the image and returns the pixel coordinates.

    M is a 3xN matrix containing the coordinates of the 3D points: M=[X;Y;Z]
    "ocam_model" contains the model of the calibrated camera.

    Returns a 2xN matrix containing the rows and columns of the points after being
    reprojected onto the image.

    Copyright (C) 2006 DAVIDE SCARAMUZZA
    Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
    """
    ss = ocam_model.ss
    xc = ocam_model.xc
    yc = ocam_model.yc
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    x, y = omni3d2pixel(ss, M)
    m = np.empty((2, len(x)))
    m[0] = x * c + y * d + xc
    m[1] = x * e + y + yc
    return m


def world2cam_fast(M, ocam_model: OCamModel):
    """Projects a 3D point on to the image and returns the pixel coordinates.
    This function uses an approximation of the inverse polynomial
    to compute the reprojected point. Therefore it is very fast.

    M is a 3xN matrix containing the coordinates of the 3D points: M=[X;Y;Z]
    "ocam_model" contains the model of the calibrated camera.

    Returns a 2xN matrix containing the rows and columns of the points after being
    reprojected onto the image.

    Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
    Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
    """
    xc = ocam_model.xc
    yc = ocam_model.yc
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    pol = ocam_model.pol
    norm = la.norm(M, axis=0)
    norm[norm == 0] = eps
    M_norm = M / norm
    theta = np.arctan(M_norm[2, :])
    rho = np.polyval(pol, theta)
    x = M_norm[0, :] * rho
    y = M_norm[1, :] * rho
    m = np.empty((2, M.shape[1]))
    m[0, :] = x * c + y * d + xc
    m[1, :] = x * e + y + yc
    return m


def omni3d2pixel(ss: np.ndarray, xx: np.ndarray):
    """Convert 3D coordinates vector into 2D pixel coordinates"""
    n_points = xx.shape[1]
    norm = la.norm(xx[:2], axis=0)
    norm[norm == 0] = eps
    xx_normed = xx / norm
    rho = np.zeros(n_points)
    poly_coef = np.repeat(ss[::-1][None], n_points, axis=0)
    poly_coef[:, -2] -= xx_normed[2]
    for j in range(n_points):
        # rho_tmp = np.roots(poly_coef[j]) # is 10x slower
        # res = np.real(rho_tmp[np.isreal(rho_tmp) & (rho_tmp > 0)])
        rho_tmp = cv2.solvePoly(poly_coef[j][::-1], maxIters=50)[1][:, 0]
        res = rho_tmp[(rho_tmp[:, 0] > 0) & (np.abs(rho_tmp[:, 1]) < 1e-10), 0]
        if len(res) == 0:
            rho[j] = np.nan
        elif len(res) > 1:
            rho[j] = np.min(res)
        else:
            rho[j] = res

    x = xx_normed[0] * rho
    y = xx_normed[1] * rho
    return x, y


def cam2world(m: np.ndarray, ocam_model: OCamModel):
    """Project a give pixel point onto the unit sphere.

    Returns the 3D coordinates of the vector emanating from the single
    effective viewpoint on the unit sphere

    Input m=[rows;cols] is a 2xN matrix containing the pixel coordinates
    of the image points.

    "ocam_model" contains the model of the calibrated camera.

    Returned M=[X;Y;Z] is a 3xN matrix with the coordinates on the unit sphere:
    thus, X^2 + Y^2 + Z^2 = 1

    Last update May 2009
    Copyright (C) 2006 DAVIDE SCARAMUZZA
    Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
    """
    ss = ocam_model.ss
    xc = ocam_model.xc
    yc = ocam_model.yc
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    A = np.array([[c, d], [e, 1]])
    T = np.array([[xc], [yc]])
    m = la.inv(A) @ (m - T)
    M = np.r_[
        m[:2, :],
        np.polyval(ss[::-1], la.norm(m[:2, :], axis=0))
    ]
    M /= la.norm(M, axis=0)  # normalizes coordinates so that they have unit length (projection onto the unit sphere)
    return M


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
    # Copyright (C) 2006 DAVIDE SCARAMUZZA
    # Author: Davide Scaramuzza - email: davsca@tiscali.it
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
