# IMUNWRAP Unwrap an omnidirectional image into a cylindrical panorama
# (cartesian to polar transformation).
#   uI = IMUNWRAP(I, center, Rmax, Rmin, bilinear)
#   I = input color image
#   center = coordinates of the circle center, center=[col;row]
#   Rmax = external radius of the omnidirectional image
#   Rmin = inner radius of the omnidirectional image

#   Copyright (C) 2006 DAVIDE SCARAMUZZA, ETH Zurich
#   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
from .libsmop import *


@function
def imunwrap(I=None, center=None, Rmax=None, Rmin=None, bilinear=None, width=None):
    varargin = imunwrap.varargin
    nargin = imunwrap.nargin

    if nargin == 3:
        Rmin = round(0.2 * Rmax)
        bilinear = 1
        width = 360
    else:
        if nargin == 4:
            bilinear = 1
            width = 360

    Rmax = round(Rmax)
    Rmin = round(Rmin)
    I = double(I)
    xc = center[2]
    yc = center[1]
    # c=round(2*np.pi*Rmax);
    c = round(width)
    cols = c
    rows = Rmax - Rmin
    Iunwraped = zeros(Rmax - Rmin, c, 3)
    R = zeros(Rmax - Rmin, c)
    G = zeros(Rmax - Rmin, c)
    B = zeros(Rmax - Rmin, c)
    RI = I[:, :, 1]
    GI = I[:, :, 2]
    BI = I[:, :, 3]
    U = zeros(size(Iunwraped))
    V = zeros(size(Iunwraped))
    r = Rmax
    J, II = meshgrid(arange(1, cols), arange(1, rows), nargout=2)
    THETA = dot(dot(- (J - 1) / width, 2), np.pi)
    RHO = Rmax + 1 - II
    X = xc + multiply(RHO, cos(THETA))
    Y = yc + multiply(RHO, sin(THETA))
    U = floor(Y)
    V = floor(X)
    IND = dot((U - 1), size(I, 1)) + V
    IND[find(logical_not((U > logical_and(1, U) < logical_and(size(I, 2), V) > logical_and(1, V) < size(I, 1))))] = 0
    fIND = find(IND)
    if bilinear:
        # if BILINEAR
        dX = Y - U
        dY = X - V
        A1 = multiply((1 - dY), (1 - dX))
        A2 = multiply((1 - dY), dX)
        A3 = multiply(dY, (1 - dX))
        A4 = multiply(dY, dX)
        indA1 = IND
        indA2 = dot((U), size(I, 1)) + V
        indA3 = dot((U - 1), size(I, 1)) + (V + 1)
        indA4 = dot(U, size(I, 1)) + (V + 1)
        R[fIND] = multiply(RI(indA1(fIND)), A1[fIND]) + \
                  multiply(RI(indA2(fIND)), A2[fIND]) + \
                  multiply(RI(indA3(fIND)), A3[fIND]) + \
                  multiply(RI(indA4(fIND)), A4[fIND])
        G[fIND] = multiply(GI(indA1(fIND)), A1[fIND]) + \
                  multiply(GI(indA2(fIND)), A2[fIND]) + \
                  multiply(GI(indA3(fIND)), A3[fIND]) + \
                  multiply(GI(indA4(fIND)), A4[fIND])
        B[fIND] = multiply(BI(indA1(fIND)), A1[fIND]) + \
                  multiply(BI(indA2(fIND)), A2[fIND]) + \
                  multiply(BI(indA3(fIND)), A3[fIND]) + \
                  multiply(BI(indA4(fIND)), A4[fIND])
    else:
        # if NOT BILINEAR
        R[find(IND)] = RI(IND(find(IND)))
        G[find(IND)] = GI(IND(find(IND)))
        B[find(IND)] = BI(IND(find(IND)))
        Iunwraped[:, :, 1] = R
        Iunwraped[:, :, 2] = G
        Iunwraped[:, :, 3] = B

    Iunwraped[:, :, 1] = R
    Iunwraped[:, :, 2] = G
    Iunwraped[:, :, 3] = B
    Iunwraped = uint8(Iunwraped)
    return Iunwraped
