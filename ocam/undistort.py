# UNDISTORT unwrap part of the image onto a plane perpendicular to the
# camera axis
#   B = UNDISTORT(OCAM_MODEL, A, FC, DISPLAY)
#   A is the input image
#   FC is a factor proportional to the distance of the camera to the plane;
#   start with FC=5 and then tune the parameter to change the result.
#   DISPLAY visualizes the output image if set to 1; its default value is
#   0.
#   B is the final image
#   Note, this function uses nearest neighbour interpolation to unwrap the
#   image point. Better undistortion methods can be implemented using
#   bilinear or bicub interpolation.
#   Note, if you want to change the size of the final image, change Nwidth
#   and Nheight
#   Author: Davide Scaramuzza, 2009

from .findinvpoly import findinvpoly
from .get_color_from_imagepoints import get_color_from_imagepoints
from .libsmop import *
from .world2cam_fast import world2cam_fast


@function
def undistort(ocam_model=None, img=None, fc=None, display=None):
    nargin = undistort.nargin

    # Parameters of the new image
    Nwidth = 640
    Nheight = 480
    Nxc = Nheight / 2
    Nyc = Nwidth / 2
    Nz = - Nwidth / fc
    if logical_not(isfield(ocam_model, 'pol')):
        width = ocam_model.width
        height = ocam_model.height
        ocam_model.pol = copy(findinvpoly(ocam_model.ss, sqrt((width / 2)**2 + (height / 2)**2)))

    if nargin < 3:
        fc = 5
        display = 0

    if length(size(img)) == 3:
        Nimg = zeros(Nheight, Nwidth, 3)
    else:
        Nimg = zeros(Nheight, Nwidth)

    i, j = meshgrid(arange(1, Nheight), arange(1, Nwidth), nargout=2)
    Nx = i - Nxc
    Ny = j - Nyc
    Nz = dot(ones(size(Nx)), Nz)
    M = concat([[ravel(Nx).T], [ravel(Ny).T], [ravel(Nz).T]])
    m = world2cam_fast(M, ocam_model)
    if length(size(img)) == 2:
        I[:, :, 1] = img
        I[:, :, 2] = img
        I[:, :, 3] = img

    r, g, b = get_color_from_imagepoints(I, m.T, nargout=3)
    Nimg = reshape(r, Nwidth, Nheight).T
    # Nimg = uint8(Nimg);
    if display:
        figure
        imagesc(Nimg)
        colormap(gray)

    # M = cam2world( distorted_points' , ocam_model );
    # M = M./(ones(3,1)*M[3,:])*(Nz);
    # 
    # ti = M[1,:] + Nxc;
    # tj = M[2,:] + Nyc;
    # 
    # scale_factor = abs(Nz);
    # 
    # und_points = [ti ; tj]';

    return Nimg
