import numpy as np
from numpy.linalg import norm, svd
from scipy.signal import convolve2d


def cornerfinder(xt: np.array, I: np.array, wintx=5, winty=5, wx2=-1, wy2=-1, line_feat=True):
    """Finds the sub-pixel corners on the image I with initial guess xt.

    xt and xc are 2xN matrices. The first component is the x coordinate
    (horizontal) and the second component is the y coordinate (vertical)
    
    Set line_feat to True to allow for extraction of line features.

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

    I = np.array(I)
    xt = np.array(xt)
    xt = np.fliplr(xt.T)

    # mask = ones(2*wintx+1,2*winty+1);
    tmp = np.exp(-(np.arange(-wintx, wintx + 1) / wintx)**2)[None]
    mask = tmp.T @ tmp
    # another mask:
    # X, Y = np.meshgrid(arange(-winty, winty), arange(-wintx, wintx))
    # mask2 = X**2 + Y**2
    # mask2[wintx + 1, winty + 1] = 1
    # mask2 = 1.0 / mask2
    # mask - mask2;

    if wx2 > 0 and wy2 > 0:
        if (wintx - wx2) >= 2 and (winty - wy2) >= 2:
            mask[wintx - wx2: wintx + 1 + wx2, winty - wy2:winty + 1 + wy2] = \
                np.zeros((2 * wx2 + 1, 2 * wy2 + 1))

    offy, offx = np.meshgrid(np.arange(-winty, winty + 1), np.arange(-wintx, wintx + 1))
    resolution = 0.005
    MaxIter = 10
    nx, ny = I.shape
    N = xt.shape[0]
    xc = xt.copy()  # first guess... they don't move !!!

    type_ = np.zeros(N)
    for i in range(N):
        v_extra = resolution + 1  # just larger than resolution
        compt = 0  # no iteration yet
        while norm(v_extra) > resolution and compt < MaxIter:
            # Coords. of the point on the initial image
            cIx = xc[i, 0]
            cIy = xc[i, 1]
            crIx = int(np.round(cIx))
            crIy = int(np.round(cIy))
            # Coefficients to compute the sub-pixel accuracy
            itIx = cIx - crIx
            itIy = cIy - crIy
            if itIx > 0:
                vIx = np.array([[itIx, 1 - itIx, 0]]).T
            else:
                vIx = np.array([[0, 1 + itIx, - itIx]]).T
            if itIy > 0:
                vIy = np.array([[itIy, 1 - itIy, 0]])
            else:
                vIy = np.array([[0, 1 + itIy, - itIy]])

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

            SI = I[xmin - 1:xmax, ymin - 1:ymax]  # The necessary neighborhood
            SI = convolve2d(convolve2d(SI, vIx, 'same'), vIy, 'same')
            SI = SI[1:2 * wintx + 4, 1:2 * winty + 4]  # The subpixel interpolated neighborhood
            gx, gy = np.gradient(np.array(SI))  # The gradient image
            gx = gx[1:2 * wintx + 2, 1:2 * winty + 2]  # extraction of the useful parts only of the gradients
            gy = gy[1:2 * wintx + 2, 1:2 * winty + 2]
            px = cIx + offx
            py = cIy + offy
            gxx = gx * gx * mask
            gyy = gy * gy * mask
            gxy = gx * gy * mask
            bb = np.array([
                (gxx * px).sum() + (gxy * py).sum(),
                (gxy * px).sum() + (gyy * py).sum()
            ])
            a = gxx.sum()
            b = gxy.sum()
            c = gyy.sum()
            dt = a * c - b**2
            xc2 = np.array([[c * bb[0] - b * bb[1], a * bb[1] - b * bb[0]]]).T / dt
            if line_feat:
                G = np.array([[a, b], [b, c]])
                U, S, Vh = svd(G)
                # If non-invertible, then project the point onto the edge orthogonal:
                if S[0] / S[1] > 50:
                    # projection operation
                    xc2 += ((xc[i, :] - xc2) * Vh[1, :]).sum(axis=0) @ Vh[1, :]
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
    delta_x = xc[:, 0] - xt[:, 0]
    delta_y = xc[:, 1] - xt[:, 1]

    bad = (np.abs(delta_x) > wintx) | (np.abs(delta_y) > winty)
    good = ~bad
    # For the diverged points, keep the original guesses:
    if np.any(bad):
        xc[bad, :] = xt[bad, :]
    xc = np.fliplr(xc)

    return xc.T, good, bad, type_
