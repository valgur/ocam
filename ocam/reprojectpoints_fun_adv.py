from .libsmop import *

from .omni3d2pixel import omni3d2pixel


@function
def reprojectpoints_fun_adv(x=None, xc=None, yc=None, ss=None, RRfin=None, ima_proc=None, Xp_abs=None, Yp_abs=None,
                            M=None, width=None, height=None):
    a = x[1]
    b = x[2]
    c = x[3]
    d = x[4]
    e = x[5]
    ssc = x(arange(6, end()))
    M[:, 3] = 1
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    for i in ima_proc.flat:
        Mc = concat([Mc, dot(RRfin[:, :, i], M.T)])
        Xpp = concat([[Xpp], [Xp_abs[:, :, i]]])
        Ypp = concat([[Ypp], [Yp_abs[:, :, i]]])

    xp, yp = omni3d2pixel(multiply(ss, ssc.T), Mc, width, height, nargout=2)
    xp = dot(xp, c) + dot(yp, d) + dot(xc, a)
    yp = dot(xp, e) + yp + dot(yc, b)
    err = sum((Xpp - xp.T)**2 + (Ypp - yp.T)**2)
    return err
