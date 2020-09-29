# WORLD2CAM projects a 3D point on to the image
#   m=WORLD2CAM(M, ocam_model) projects a 3D point on to the
#   image and returns the pixel coordinates.
#   
#   M is a 3xN matrix containing the coordinates of the 3D points: M=[X;Y;Z]
#   "ocam_model" contains the model of the calibrated camera.
#   m=[rows;cols] is a 2xN matrix containing the returned rows and columns of the points after being
#   reproject onto the image.
#   
#   Copyright (C) 2006 DAVIDE SCARAMUZZA   
#   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org


from .calib_data import OCamModel
from .libsmop import *
from .omni3d2pixel import omni3d2pixel


@function
def world2cam(M, ocam_model: OCamModel):
    ss = ocam_model.ss
    xc = ocam_model.xc
    yc = ocam_model.yc
    width = ocam_model.width
    height = ocam_model.height
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    x, y = omni3d2pixel(ss, M, width, height, nargout=2)
    m[1, :] = dot(x, c) + dot(y, d) + xc
    m[2, :] = dot(x, e) + y + yc
    return m
