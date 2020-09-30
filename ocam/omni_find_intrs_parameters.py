from .libsmop import *


@function
def omni_find_intrs_parameters(taylor_order=None, xc=None, yc=None, ima_proc=None, Xp_abs=None, Yp_abs=None, Xt=None,
                               Yt=None, RRfin=None):
    varargin = omni_find_intrs_parameters.varargin
    nargin = omni_find_intrs_parameters.nargin

    # Find the other parameters
    Xp = Xp_abs - xc
    Yp = Yp_abs - yc
    PP = copy([])
    QQ = copy([])
    for i in ima_proc.flat:
        RRdef = RRfin[:, :, i]
        R11 = RRdef[1, 1]
        R21 = RRdef[2, 1]
        R31 = RRdef[3, 1]
        R12 = RRdef[1, 2]
        R22 = RRdef[2, 2]
        R32 = RRdef[3, 2]
        T1 = RRdef[1, 3]
        T2 = RRdef[2, 3]
        T3 = RRdef[3, 3]
        Xpt = Xp[:, :, i]
        Ypt = Yp[:, :, i]
        MA = multiply(R21, Xt) + multiply(R22, Yt) + T2
        MB = multiply(Ypt, (multiply(R31, Xt) + multiply(R32, Yt) + T3))
        MC = multiply(R11, Xt) + multiply(R12, Yt) + T1
        MD = multiply(Xpt, (multiply(R31, Xt) + multiply(R32, Yt) + T3))
        rho = copy([])
        for j in arange(2, taylor_order).flat:
            rho[:, :, j] = (sqrt(Xpt**2 + Ypt**2))**j
        PP1 = concat([[MA], [MC]])
        for j in arange(2, taylor_order).flat:
            PP1 = concat(
                [PP1, concat([[multiply(MA, rho[:, :, j])], [multiply(MC, rho[:, :, j])]])])
        PP = concat([[PP], [PP1]])
        QQ = concat([[QQ], [MB], [MD]])

    ss = dot(pinv(PP), QQ)
    ss = concat([[ss[1]], [0], [ss(arange(2, end()))]])

    return ss
