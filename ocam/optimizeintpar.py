from .calib_data import CalibData
from .libsmop import *
from .prova2 import prova2
from .reprojectpoints_adv import reprojectpoints_adv


@function
def optimizeintpar(calib_data: CalibData):
    options = optimset('Display', 'final', 'LargeScale', 'on', 'TolX', 0.0001, 'TolFun', 0.0001, 'DerivativeCheck',
                       'off', 'Diagnostics', 'off', 'Jacobian', 'off', 'JacobMult', [], 'JacobPattern',
                       'sparse(ones(Jrows,Jcols))', 'MaxFunEvals', '100*numberOfVariables', 'DiffMaxChange', 0.1,
                       'DiffMinChange', 1e-08, 'PrecondBandWidth', 0, 'TypicalX', 'ones(numberOfVariables,1)',
                       'MaxPCGIter', 'max(1,floor(numberOfVariables/2))', 'TolPCG', 0.1, 'MaxIter', 10000,
                       'LineSearchType', 'quadcubic', 'LevenbergMarquardt', 'on')
    if (logical_and(logical_and(isempty(calib_data.ocam_model.c), isempty(calib_data.ocam_model.d)),
                    isempty(calib_data.ocam_model.e))):
        calib_data.ocam_model.c = 1
        calib_data.ocam_model.d = 0
        calib_data.ocam_model.e = 0

    int_par = concat(
        [calib_data.ocam_model.c, calib_data.ocam_model.d, calib_data.ocam_model.e, calib_data.ocam_model.xc,
         calib_data.ocam_model.yc])
    ssout, resnorm, residual, exitflag, output = lsqnonlin(prova2, int_par, - inf, inf, options,
                                                           calib_data.ocam_model.ss, calib_data.RRfin,
                                                           calib_data.ima_proc, calib_data.Xp_abs, calib_data.Yp_abs, M,
                                                           calib_data.ocam_model.width, calib_data.ocam_model.height,
                                                           nargout=5)
    print('Camera model optimized')
    c = ssout[1]
    d = ssout[2]
    e = ssout[3]
    xc = ssout[4]
    yc = ssout[5]
    # calib_data.ocam_model.ss=ss;
    calib_data.ocam_model.xc = xc
    calib_data.ocam_model.yc = yc
    calib_data.ocam_model.c = c
    calib_data.ocam_model.d = d
    calib_data.ocam_model.e = e
    # calib_data.ocam_model.width=width;
    # calib_data.ocam_model.height=height;
    reprojectpoints_adv(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc, calib_data.Xp_abs,
                        calib_data.Yp_abs, M)
    print(calib_data.ocam_model.ss)
