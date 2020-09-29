# FINDINVPOLY finds the inverse polynomial specified in the argument.
#   [POL, ERR, N] = FINDINVPOLY(SS, RADIUS, N) finds an approximation of the inverse polynomial specified in OCAM_MODEL.SS.
#   The returned polynomial POL is used in WORLD2CAM_FAST to compute the reprojected point very efficiently.
#   
#   SS is the polynomial which describe the mirrror/lens model.
#   RADIUS is the radius (pixels) of the omnidirectional picture.
#   ERR is the error (pixel) that you commit in using the returned
#   polynomial instead of the inverse SS. N is searched so that
#   that ERR is < 0.01 pixels.

#   Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
#   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org


from .libsmop import *


@function
def findinvpolyUrban(ss=None, I=None):
    nargin = findinvpolyUrban.nargin

    disp('choose about 10 points along the border of the visible image region! Press enter if finished')
    figure
    imshow(uint8(I))
    x, y = ginput
    X = concat([x.T, y.T])
    _, A, B, _ = fitellipse(X, nargout=4)
    radius = mean(concat([A, B]))
    if nargin < 3:
        maxerr = inf
        N = 1
        while maxerr > 0.01:
            N += 1
            pol, err, N = findinvpoly2(ss, radius, N, nargout=3)
            maxerr = max(err)

    else:
        pol, err, N = findinvpoly2(ss, radius, N, nargout=3)
    return pol, err, N


@function
def findinvpoly2(ss=None, radius=None, N=None):
    theta = arange(- np.pi / 2, 1.4, 0.01)
    r = invFUN(ss, theta, radius)
    ind = find(r != inf)
    theta = theta(ind)
    r = r(ind)
    pol = polyfit(theta, r, N)
    err = abs(r - polyval(pol, theta))
    return pol, err, N


@function
def invFUN(ss=None, theta=None, radius=None):
    m = np.tan(theta)
    r = copy([])
    poly_coef = ss(arange(end(), 1, - 1))
    poly_coef_tmp = copy(poly_coef)
    for j in arange(1, length(m)).flat:
        poly_coef_tmp[end() - 1] = poly_coef(end() - 1) - m(j)
        rhoTmp = roots(poly_coef_tmp)
        res = rhoTmp(find(imag(rhoTmp) == logical_and(0, rhoTmp) > logical_and(0, rhoTmp) < radius))
        if logical_or(isempty(res), length(res)) > 1:
            r[j] = inf
        else:
            r[j] = res
    return r
