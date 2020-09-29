from .libsmop import *

from .omni3d2pixel import omni3d2pixel


@function
def reprojectpoints_fun(Xt=None, Yt=None, Xp_abs=None, Yp_abs=None, xc=None, yc=None, RRfin=None, ss=None,
                        ima_proc=None, width=None, height=None):
    err = copy([])
    stderr = copy([])
    MSE = 0
    for i in ima_proc.flat:
        xx = dot(RRfin[:, :, i], concat([[Xt.T], [Yt.T], [ones(size(Xt.T))]]))
        Xp_reprojected, Yp_reprojected = omni3d2pixel(ss, xx, width, height)
        if np.isnan(Xp_reprojected):
            return np.nan
        stt = sqrt(
            (Xp_abs[:, :, i] - xc - Xp_reprojected.T)**2 +
            (Yp_abs[:, :, i] - yc - Yp_reprojected.T)**2
        )
        err[i] = np.mean(stt)
        stderr[i] = np.std(stt)
        MSE += sum(
            (Xp_abs[:, :, i] - xc - Xp_reprojected.T)**2 +
            (Yp_abs[:, :, i] - yc - Yp_reprojected.T)**2
        )
    return MSE
