import numpy as np
from numpy import linalg as la

from .calib_data import OCamModel


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
