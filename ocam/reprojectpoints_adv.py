from .libsmop import *

from .world2cam import world2cam


@function
def reprojectpoints_adv(ocam_model=None, RRfin=None, ima_proc=None, Xp_abs=None, Yp_abs=None, M=None):
    width = ocam_model.width
    height = ocam_model.height
    c = ocam_model.c
    d = ocam_model.d
    e = ocam_model.e
    ss = ocam_model.ss
    xc = ocam_model.xc
    yc = ocam_model.yc
    M[:, 3] = 1
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    count = 0
    MSE = 0
    for i in ima_proc.flat:
        count += 1
        Mc = dot(RRfin[:, :, i], M.T)
        #     [xp,yp]=omni3d2pixel(ss,Mc, width, height);
        #     xp=xp*c + yp*d + xc;
        #     yp=xp*e + yp + yc;
        m = world2cam(Mc, ocam_model)
        xp = m[1, :]
        yp = m[2, :]
        sqerr = (Xp_abs[:, :, i] - xp.T)**2 + (Yp_abs[:, :, i] - yp.T)**2
        err[count] = mean(sqrt(sqerr))
        stderr[count] = std(sqrt(sqerr))
        MSE += sum(sqerr)

    print('\n Average reprojection error computed for each chessboard [pixels]:\n')
    for i in arange(1, length(err)).flat:
        fprintf(' %3.2f %c %3.2f\n', err[i], 177, stderr[i])

    print('\n Average error [pixels]\n\n %f\n', mean(err), end='')
    print('\n Sum of squared errors\n\n %f\n', MSE, end='')

    return err, stderr, MSE
