from .libsmop import *

from .omni3d2pixel import omni3d2pixel


@function
def reprojectpoints_fun_adv_all(x=None, xc=None, yc=None, ss=None, RRfin=None, ima_proc=None, Xp_abs=None, Yp_abs=None,
                                M=None, width=None, height=None):
    c = x[1]
    d = x[2]
    e = x[3]
    ssc = x(arange(end() - length(ss) + 1, end()))
    M[:, 3] = 1
    Mc = copy([])
    Xpp = copy([])
    Ypp = copy([])
    count = 0
    MSE = 0
    for i in ima_proc.flat:
        count += 1
        Mc = dot(RRfin[:, :, i], M.T)
        xp, yp = omni3d2pixel(multiply(ss, ssc), Mc, width, height, nargout=2)
        xp = dot(xp, c) + dot(yp, d) + xc
        yp = dot(xp, e) + yp + yc
        sqerr = (Xp_abs[:, :, i] - xp.T)**2 + (Yp_abs[:, :, i] - yp.T)**2
        err[count] = mean(sqrt(sqerr))
        stderr[count] = std(sqrt(sqerr))
        MSE += sum(sqerr)

    print('\n Average reprojection error computed for each chessboard [pixels]:\n')
    for i in arange(1, length(err)).flat:
        fprintf(' %3.2f %c %3.2f\n', err(i), 177, stderr(i))

    print('\n Average error [pixels]\n\n %f\n', mean(err), end='')
    print('\n Sum of squared errors\n\n %f\n', MSE, end='')
    return err
