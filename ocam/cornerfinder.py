from numpy import meshgrid, fliplr
from numpy.linalg import norm

from .libsmop import *


@function
def cornerfinder(xt=None, I=None, wintx=None, winty=None, wx2=None, wy2=None):
    # [xc] = cornerfinder(xt,I);
    #
    # Finds the sub-pixel corners on the image I with initial guess xt
    # xt and xc are 2xN matrices. The first component is the x coordinate
    # (horizontal) and the second component is the y coordinate (vertical)
    # 
    # Based on Harris corner finder method
    #
    # Finds corners to a precision below .1 pixel!
    # Oct. 14th, 1997 - UPDATED to work with vertical and horizontal edges as well!!!
    # Sept 1998 - UPDATED to handle diverged points: we keep the original points
    # good is a binary vector indicating wether a feature point has been properly
    # found.
    #
    # Add a zero zone of size wx2,wy2
    # July 15th, 1999 - Bug on the mask building... fixed + change to Gaussian mask with higher
    # resolution and larger number of iterations.
    #
    # California Institute of Technology
    # (c) Jean-Yves Bouguet -- Oct. 14th, 1997

    varargin = cornerfinder.varargin
    nargin = cornerfinder.nargin

    line_feat = 1

    xt = xt.T
    xt = fliplr(xt)
    if nargin < 4:
        winty = 5
        if nargin < 3:
            wintx = 5

    if nargin < 6:
        wx2 = - 1
        wy2 = - 1

    # mask = ones(2*wintx+1,2*winty+1);
    mask = exp(- ((arange(- wintx, wintx)).T / wintx)**2) @ exp(- ((arange(- winty, winty)) / winty)**2)
    # another mask:
    X, Y = meshgrid(arange(- winty, winty), arange(- wintx, wintx))
    mask2 = X**2 + Y**2
    mask2[wintx + 1, winty + 1] = 1
    mask2 = 1.0 / mask2
    # mask - mask2;

    if logical_and((wx2 > 0), (wy2 > 0)):
        if logical_and(((wintx - wx2) >= 2), ((winty - wy2) >= 2)):
            mask[arange(wintx + 1 - wx2, wintx + 1 + wx2), arange(winty + 1 - wy2, winty + 1 + wy2)] = zeros(
                dot(2, wx2) + 1, dot(2, wy2) + 1)

    offx = dot(concat([arange(- wintx, wintx)]).T, ones(1, dot(2, winty) + 1))
    offy = dot(ones(dot(2, wintx) + 1, 1), concat([arange(- winty, winty)]))
    resolution = 0.005
    MaxIter = 10
    nx, ny = I.shape
    N = size(xt, 1)
    xc = xt

    type_ = zeros(1, N)
    for i in arange(1, N).flat:
        v_extra = resolution + 1
        compt = 0
        while (norm(v_extra) > resolution) and (compt < MaxIter):
            cIx = xc(i, 1)
            cIy = xc(i, 2)
            crIx = round(cIx)
            crIy = round(cIy)
            itIx = cIx - crIx
            itIy = cIy - crIy
            if itIx > 0:
                vIx = concat([itIx, 1 - itIx, 0]).T
            else:
                vIx = concat([0, 1 + itIx, - itIx]).T
            if itIy > 0:
                vIy = concat([itIy, 1 - itIy, 0])
            else:
                vIy = concat([0, 1 + itIy, - itIy])
            if crIx - wintx - 2 < 1:
                xmin = 1
                xmax = dot(2, wintx) + 5
            else:
                if crIx + wintx + 2 > nx:
                    xmax = nx
                    xmin = nx - dot(2, wintx) - 4
                else:
                    xmin = crIx - wintx - 2
                    xmax = crIx + wintx + 2
            if crIy - winty - 2 < 1:
                ymin = 1
                ymax = dot(2, winty) + 5
            else:
                if crIy + winty + 2 > ny:
                    ymax = ny
                    ymin = ny - dot(2, winty) - 4
                else:
                    ymin = crIy - winty - 2
                    ymax = crIy + winty + 2
            SI = I(arange(xmin, xmax), arange(ymin, ymax))
            SI = conv2(conv2(SI, vIx, 'same'), vIy, 'same')
            SI = SI(arange(2, dot(2, wintx) + 4), arange(2, dot(2, winty) + 4))
            gy, gx = gradient(SI, nargout=2)
            gx = gx(arange(2, dot(2, wintx) + 2), arange(2, dot(2, winty) + 2))
            gy = gy(arange(2, dot(2, wintx) + 2), arange(2, dot(2, winty) + 2))
            px = cIx + offx
            py = cIy + offy
            gxx = multiply(multiply(gx, gx), mask)
            gyy = multiply(multiply(gy, gy), mask)
            gxy = multiply(multiply(gx, gy), mask)
            bb = concat(
                [[sum(sum(multiply(gxx, px) + multiply(gxy, py)))], [sum(sum(multiply(gxy, px) + multiply(gyy, py)))]])
            a = sum(sum(gxx))
            b = sum(sum(gxy))
            c = sum(sum(gyy))
            dt = dot(a, c) - b**2
            xc2 = concat([dot(c, bb[1]) - dot(b, bb[2]), dot(a, bb[2]) - dot(b, bb[1])]) / dt
            if line_feat:
                G = concat([[a, b], [b, c]])
                U, S, V = svd(G, nargout=3)
                # If non-invertible, then project the point onto the edge orthogonal:
                if S(1, 1) / S(2, 2) > 50:
                    xc2 += dot(sum(multiply((xc[i, :] - xc2), V[:, 2].T)), V[:, 2].T)
                    type_[i] = 1
            #      G = [a b;b c];
            #      [U,S,V]  = svd(G);
            #      if S(1,1)/S(2,2) > 150,
            #	 bb2 = U'*bb;
            #	 xc2 = (V*[bb2[1]/S(1,1) ;0])';
            #      else
            #	 xc2 = [c*bb[1]-b*bb[2] a*bb[2]-b*bb[1]]/dt;
            #      end;
            # if (abs(a)> 50*abs(c)),
            #	 xc2 = [(c*bb[1]-b*bb[2])/dt xc(i,2)];
            #      elseif (abs(c)> 50*abs(a))
            #	 xc2 = [xc(i,1) (a*bb[2]-b*bb[1])/dt];
            #      else
            #	 xc2 = [c*bb[1]-b*bb[2] a*bb[2]-b*bb[1]]/dt;
            #      end;
            # keyboard;
            v_extra = xc[i, :] - xc2
            xc[i, :] = xc2
            #      keyboard;
            compt += 1

    # check for points that diverge:

    delta_x = xc[:, 1] - xt[:, 1]
    delta_y = xc[:, 2] - xt[:, 2]
    # keyboard;

    bad = logical_or((abs(delta_x) > wintx), (abs(delta_y) > winty))
    good = logical_not(bad)
    in_bad = find(bad)
    # For the diverged points, keep the original guesses:

    xc[in_bad, :] = xt[in_bad, :]
    xc = fliplr(xc)
    xc = xc.T
    bad = bad.T
    good = good.T

    return xc, good, bad, type
