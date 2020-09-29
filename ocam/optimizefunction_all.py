from .calib_data import CalibData
from .libsmop import *
from .reprojectpoints_adv import reprojectpoints_adv
from .rodrigues import rodrigues


@function
def optimizefunction_all(calib_data: CalibData):
    print('\nThis function refines all calibration parameters (both EXTRINSIC and INTRINSIC)')
    print('by using a non linear minimization method ')
    print('Because of the computations involved this refinement can take several minutes')
    print('\nWARNING: Search space dimension increases linearly with number of calibration images\n')
    print('Press Enter to continue or Ctrl+C to abort now')
    pause
    options = optimset('Display', 'iter', 'LargeScale', 'off', 'TolX', 0.0001, 'TolFun', 0.0001, 'DerivativeCheck',
                       'off', 'Diagnostics', 'off', 'Jacobian', 'off', 'JacobMult', [], 'JacobPattern',
                       'sparse(ones(Jrows,Jcols))', 'MaxFunEvals', '100*numberOfVariables', 'DiffMaxChange', 0.1,
                       'DiffMinChange', 1e-08, 'PrecondBandWidth', 0, 'TypicalX', 'ones(numberOfVariables,1)',
                       'MaxPCGIter', 'max(1,floor(numberOfVariables/2))', 'TolPCG', 0.1, 'MaxIter', 10000, 'Algorithm',
                       'trust-region-reflective')
    M = concat([calib_data.Xt, calib_data.Yt, zeros(size(calib_data.Xt))])

    # costruisci vettore di stato
    if logical_and(
            logical_and(logical_not(isempty(calib_data.ocam_model.c)), logical_not(isempty(calib_data.ocam_model.d))),
            logical_not(isempty(calib_data.ocam_model.e))):
        x0 = concat([[calib_data.ocam_model.c], [calib_data.ocam_model.d], [calib_data.ocam_model.e]])
    else:
        x0 = concat([[1], [1], [1]])

    for i in calib_data.ima_proc.flat:
        r1 = calib_data.RRfin[:, 1, i]
        r2 = calib_data.RRfin[:, 2, i]
        rod = rodrigues(concat([r1, r2, cross(r1, r2)]))
        Tod = calib_data.RRfin[:, 3, i]
        x0 = concat([[x0], [rod], [Tod]])

    ss0 = calib_data.ocam_model.ss
    x0 = concat([[x0], [calib_data.ocam_model.xc], [calib_data.ocam_model.yc]])
    x0 = concat([[x0], [ones(size(ss0))]])
    tic
    allout, resnorm, residual, exitflag, output = lsqnonlin(prova_all, x0, [], [], options, ss0, calib_data.ima_proc,
                                                            calib_data.Xp_abs, calib_data.Yp_abs, M,
                                                            calib_data.ocam_model.width, calib_data.ocam_model.height,
                                                            nargout=5)
    toc
    xc = allout(end() - length(calib_data.ocam_model.ss) - 1)
    yc = allout(end() - length(calib_data.ocam_model.ss))
    c = allout[1]
    d = allout[2]
    e = allout[3]
    # costruisci vettore ssc
    ssc = allout(arange(end() - length(ss0) + 1, end()))
    ss = multiply(ss0, ssc)
    # costruisci RRfin
    count = 0
    for i in calib_data.ima_proc.flat:
        Rod = rodrigues(allout(arange(dot(6, count) + 4, dot(6, count) + 6)))
        Tod = allout(arange(dot(6, count) + 7, dot(6, count) + 9))
        RRfinOpt[:, :, i] = Rod
        RRfinOpt[:, 3, i] = Tod
        count += 1

    calib_data.RRfin = copy(RRfinOpt)
    calib_data.ocam_model.ss = ss
    calib_data.ocam_model.xc = xc
    calib_data.ocam_model.yc = yc
    calib_data.ocam_model.c = c
    calib_data.ocam_model.d = d
    calib_data.ocam_model.e = e
    reprojectpoints_adv(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc, calib_data.Xp_abs,
                        calib_data.Yp_abs, M)
    ss
