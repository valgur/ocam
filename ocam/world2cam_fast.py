# WORLD2CAM projects a 3D point on to the image
#   m=WORLD2CAM_FAST(M, ocam_model) projects a 3D point on to the
#   image and returns the pixel coordinates. This function uses an approximation of the inverse
#   polynomial to compute the reprojected point. Therefore it is very fast.
#   
#   M is a 3xN matrix containing the coordinates of the 3D points: M=[X;Y;Z]
#   "ocam_model" contains the model of the calibrated camera.
#   m=[rows;cols] is a 2xN matrix containing the returned rows and columns of the points after being
#   reproject onto the image.
#   
#   Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich  
#   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org


from .calib_data import OCamModel
from .libsmop import *


@function
def world2cam_fast(M, ocam_model: OCamModel):
    ss = ocam_model.ss
    xc = ocam_model.xc
    yc = ocam_model.yc
    width = ocam_model.width
    height = ocam_model.height
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    pol = ocam_model.pol
    npoints = size(M, 2)
    theta = zeros(1, npoints)
    NORM = sqrt(M[1, :]**2 + M[2, :]**2)
    ind0 = find(NORM == 0)

    NORM[ind0] = eps

    theta = atan(M[3, :] / NORM)
    rho = polyval(pol, theta)

    x = multiply(M[1, :] / NORM, rho)
    y = multiply(M[2, :] / NORM, rho)
    # Add center coordinates
    m[1, :] = dot(x, c) + dot(y, d) + xc
    m[2, :] = dot(x, e) + y + yc
    return m
