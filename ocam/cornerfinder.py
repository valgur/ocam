from numpy.linalg import norm

from .libsmop import *


def cornerfinder(xt: matlabarray, I: matlabarray, wintx=5, winty=5, wx2=-1, wy2=-1, line_feat=True):
    """Finds the sub-pixel corners on the image I with initial guess xt.

    xt and xc are 2xN matrices. The first component is the x coordinate
    (horizontal) and the second component is the y coordinate (vertical)

    Based on Harris corner finder method

    Finds corners to a precision below .1 pixel!
    Oct. 14th, 1997 - UPDATED to work with vertical and horizontal edges as well!!!
    Sept 1998 - UPDATED to handle diverged points: we keep the original points
    good is a binary vector indicating whether a feature point has been properly found.

    Add a zero zone of size wx2,wy2
    July 15th, 1999 - Bug on the mask building... fixed + change to Gaussian mask with higher
    resolution and larger number of iterations.

    California Institute of Technology
    (c) Jean-Yves Bouguet -- Oct. 14th, 1997
    """

    xt = fliplr(xt.T)

    # mask = ones(2*wintx+1,2*winty+1);
    tmp = exp(-(arange(-wintx, wintx) / wintx)**2)
    mask = tmp.T @ tmp
    # another mask:
    # X, Y = np.meshgrid(arange(-winty, winty), arange(-wintx, wintx))
    # mask2 = X**2 + Y**2
    # mask2[wintx + 1, winty + 1] = 1
    # mask2 = 1.0 / mask2
    # mask - mask2;

    if wx2 > 0 and wy2 > 0:
        if (wintx - wx2) >= 2 and (winty - wy2) >= 2:
            mask[wintx + 1 - wx2: wintx + 1 + wx2, winty + 1 - wy2:winty + 1 + wy2] = \
                zeros(2 * wx2 + 1, 2 * wy2 + 1)

    offy, offx = np.meshgrid(arange(-winty, winty), arange(-wintx, wintx))
    resolution = 0.005
    MaxIter = 10
    nx, ny = I.shape
    N = xt.shape[0]
    xc = copy(xt)

    type_ = zeros(1, N)
    for i in range(1, N + 1):
        v_extra = resolution + 1
        compt = 0
        while norm(v_extra) > resolution and compt < MaxIter:
            cIx = float(xc[i, 1])
            cIy = float(xc[i, 2])
            crIx = np.round(cIx)
            crIy = np.round(cIy)
            itIx = cIx - crIx
            itIy = cIy - crIy
            if itIx > 0:
                vIx = copy([itIx, 1 - itIx, 0]).T
            else:
                vIx = copy([0, 1 + itIx, - itIx]).T
            if itIy > 0:
                vIy = copy([itIy, 1 - itIy, 0])
            else:
                vIy = copy([0, 1 + itIy, - itIy])

            # What if the sub image is not in?
            if crIx - wintx - 2 < 1:
                xmin = 1
                xmax = 2 * wintx + 5
            elif crIx + wintx + 2 > nx:
                xmax = nx
                xmin = nx - 2 * wintx - 4
            else:
                xmin = crIx - wintx - 2
                xmax = crIx + wintx + 2
            if crIy - winty - 2 < 1:
                ymin = 1
                ymax = 2 * winty + 5
            elif crIy + winty + 2 > ny:
                ymax = ny
                ymin = ny - 2 * winty - 4
            else:
                ymin = crIy - winty - 2
                ymax = crIy + winty + 2

            SI = I[xmin:xmax, ymin:ymax]
            SI = conv2(conv2(SI, vIx, 'same'), vIy, 'same')
            SI = SI[2:2 * wintx + 4, 2:2 * winty + 4]
            gx, gy = np.gradient(np.array(SI))
            gx = matlabarray(gx)[2:2 * wintx + 2, 2:2 * winty + 2]
            gy = matlabarray(gy)[2:2 * wintx + 2, 2:2 * winty + 2]
            px = cIx + offx
            py = cIy + offy
            gxx = gx * gx * mask
            gyy = gy * gy * mask
            gxy = gx * gy * mask
            bb = concat([
                (gxx * px).sum() + (gxy * py).sum(),
                (gxy * px).sum() + (gyy * py).sum()
            ])
            a = gxx.sum()
            b = gxy.sum()
            c = gyy.sum()
            dt = a * c - b**2
            xc2 = concat([c * bb[1] - b * bb[2], a * bb[2] - b * bb[1]]).T / dt
            if line_feat:
                G = matlabarray([[a, b], [b, c]])
                U, S, V = svd(G)
                # If non-invertible, then project the point onto the edge orthogonal:
                if S[1, 1] / S[2, 2] > 50:
                    xc2 += sum(multiply((xc[i, :] - xc2), V[:, 2].T)) @ V[:, 2].T
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
                #	 xc2 = [(c*bb[1]-b*bb[2])/dt xc[i,2]];
                #      elseif (abs(c)> 50*abs(a))
                #	 xc2 = [xc[i,1] (a*bb[2]-b*bb[1])/dt];
                #      else
                #	 xc2 = [c*bb[1]-b*bb[2] a*bb[2]-b*bb[1]]/dt;
                #      end;
            v_extra = xc[i, :] - xc2
            xc[i, :] = np.squeeze(xc2)
            compt += 1

    # check for points that diverge:
    delta_x = xc[:, 1] - xt[:, 1]
    delta_y = xc[:, 2] - xt[:, 2]

    bad = logical_or(abs(delta_x) > wintx, abs(delta_y) > winty)
    good = logical_not(bad)
    in_bad = find(bad)
    # For the diverged points, keep the original guesses:
    if len(in_bad):
        xc[in_bad, :] = xt[in_bad, :]
    xc = fliplr(xc)
    xc = xc.T
    bad = bad.T
    good = good.T

    return xc, good, bad, type_
