###########################################################################
#   Copyright (C) 2006 DAVIDE SCARAMUZZA
#
#   Author: Davide Scaramuzza - email: davsca@tiscali.it
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#   USA
############################################################################
#
# 04.03.2014 by Steffen Urban
# this is a modified file from
# Davide Scaramuzzas Toolbox OcamCalib
# filename: optimizefunction.m
# Code was added to calculate the standard deviations of the ext ori
from .calib_data import CalibData
from .libsmop import *

from .reprojectpoints_adv import reprojectpoints_adv
from .reprojectpoints_quiet import reprojectpoints_quiet
from .rodrigues import rodrigues


@function
def optimizefunction(calib_data: CalibData):
    tol_MSE = 0.0001
    MSE_old = 0
    MSE_new = np.inf
    iter = 0
    print('This function alternately refines EXTRINSIC and INTRINSIC calibration parameters')
    print('by using an interative non linear minimization method ')
    print('Because of the computations involved this refinement can take some seconds')
    print('Loop interrupting: Press enter to stop refinement. (OCamCalib GUI must be selected!)')
    # max_iter = input('\n Maximum number of iterations ([] = 100, 0 = abort, -1 = no limit) = ');

    # if ~isempty(max_iter)
    #     if max_iter == 0
    #         return;
    #     elseif max_iter == -1;
    #         max_iter = np.inf;
    #     end
    # else
    max_iter = 100
    # end

    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    pause(0.01)
    figure[1]

    pause(0.01)
    set(gcf, 'CurrentChar', 'a')
    while (iter < max_iter and abs(MSE_new - MSE_old) > tol_MSE) and get(gcf, 'CurrentChar') != 13:

        iter += 1
        print('\nIteration %i\n', iter, end='')
        print('Starting refinement of EXTRINSIC parameters...')
        options = optimset('Display', 'off', 'LargeScale', 'off', 'TolX', 0.0001, 'TolFun', 0.0001, 'DerivativeCheck',
                           'off', 'Diagnostics', 'off', 'Jacobian', 'off', 'JacobMult', [], 'JacobPattern',
                           'sparse(ones(Jrows,Jcols))', 'MaxFunEvals', '100*numberOfVariables', 'DiffMaxChange', 0.1,
                           'DiffMinChange', 1e-08, 'PrecondBandWidth', 0, 'TypicalX', 'ones(numberOfVariables,1)',
                           'MaxPCGIter', 'max(1,floor(numberOfVariables/2))', 'TolPCG', 0.1, 'MaxIter', 10000)
        if (logical_and(logical_and(isempty(calib_data.ocam_model.c), isempty(calib_data.ocam_model.d)),
                        isempty(calib_data.ocam_model.e))):
            calib_data.ocam_model.c = 1
            calib_data.ocam_model.d = 0
            calib_data.ocam_model.e = 0
        int_par = concat(
            [calib_data.ocam_model.c, calib_data.ocam_model.d, calib_data.ocam_model.e, calib_data.ocam_model.xc,
             calib_data.ocam_model.yc])
        M = concat([calib_data.Xt, calib_data.Yt, zeros(size(calib_data.Xt))])
        print('Optimizing chessboard pose ', end='')
        for i in calib_data.ima_proc.flat:
            print('%d,  ', i, end='')
            R = calib_data.RRfin[:, :, i]
            R[:, 3] = cross(R[:, 1], R[:, 2])
            r = rodrigues(R)
            t = calib_data.RRfin[:, 3, i]
            x0 = concat([r[1], r[2], r[3], t[1], t[2], t[3]])
            x0, __, vExtr, __, __, __, jacExtr = lsqnonlin(prova, x0, dot(- inf, ones(size(x0))),
                                                           dot(inf, ones(size(x0))), options, calib_data.ocam_model.ss,
                                                           int_par, calib_data.Xp_abs[:, :, i],
                                                           calib_data.Yp_abs[:, :, i], M,
                                                           calib_data.ocam_model.width, calib_data.ocam_model.height,
                                                           nargout=7)
            RRfinOpt[:, :, i] = rodrigues(x0(arange(1, 3)))
            RRfinOpt[:, 3, i] = x0(arange(4, 6)).T
            #  added code, steffen urban
            sigma0q = 1**2
            rows, cols = size(jacExtr, nargout=2)
            v = copy(vExtr)
            J = copy(jacExtr)
            Qll = dot(sigma0q, eye(size(J, 1), size(J, 1)))
            P0 = inv(Qll)
            calib_data.statEO[i].sg0 = copy((dot(dot(v.T, P0), v)) / (rows - cols))
            calib_data.statEO[i].Exx = copy(dot(calib_data.statEO[i].sg0, pinv(dot(dot(J.T, P0), J))))
            calib_data.statEO[i].varEO = copy(diag(calib_data.statEO[i].Exx))
            calib_data.statEO[i].stdEO = copy(sqrt(calib_data.statEO[i].varEO))
            calib_data.optimized = true
        calib_data.RRfin = copy(RRfinOpt)
        print('\nStarting refinement of INTRINSIC parameters...')
        ss0 = calib_data.ocam_model.ss
        #     options=optimset('Display','off',...
        #         'LargeScale','off', ...
        #         'TolX',1e-4,...
        #         'TolFun',1e-4,...
        #         'DerivativeCheck','off',...
        #         'Diagnostics','off',...
        #         'Jacobian','off',...
        #         'JacobMult',[],... # JacobMult set to [] by default
        #         'JacobPattern','sparse(ones(Jrows,Jcols))',...
        #         'MaxFunEvals','100*numberOfVariables',...
        #         'DiffMaxChange',1e-1,...
        #         'DiffMinChange',1e-8,...
        #         'PrecondBandWidth',0,...
        #         'TypicalX','ones(numberOfVariables,1)',...
        #         'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
        #         'TolPCG',0.1,...
        #         'MaxIter',10000);
        options = optimset('Display', 'off', 'LargeScale', 'off', 'TolX', 0.0001, 'TolFun', 0.0001, 'Algorithm',
                           'levenberg-marquardt')
        f0 = concat([1, 1, calib_data.ocam_model.c, calib_data.ocam_model.d, calib_data.ocam_model.e,
                     ones(1, size(calib_data.ocam_model.ss, 1))])
        lb = concat([0, 0, 0, - 1, - 1, zeros(1, size(calib_data.ocam_model.ss, 1))])
        ub = concat([2, 2, 2, 1, 1, dot(2, ones(1, size(calib_data.ocam_model.ss, 1)))])
        #     [ssout,~,vIO,~,~,~,jacsIO] =lsqnonlin(@prova3,f0,lb,ub,options,calib_data.ocam_model.xc,calib_data.ocam_model.yc,ss0,calib_data.RRfin,calib_data.ima_proc,calib_data.Xp_abs,calib_data.Yp_abs,M, calib_data.ocam_model.width, calib_data.ocam_model.height);
        ssout, __, vIO, __, __, __, jacsIO = lsqnonlin(prova3, f0, dot(- inf, ones(1, length(x0))),
                                                       dot(inf, ones(1, length(x0))), options, calib_data.ocam_model.xc,
                                                       calib_data.ocam_model.yc, ss0, calib_data.RRfin,
                                                       calib_data.ima_proc, calib_data.Xp_abs, calib_data.Yp_abs, M,
                                                       calib_data.ocam_model.width, calib_data.ocam_model.height,
                                                       nargout=7)
        #  added code , steffen urban
        sigma0q = 1**2
        rows, cols = size(jacsIO, nargout=2)
        jacsIO[:, 7] = copy([])
        v = copy(vIO)
        J = copy(jacsIO)
        Qll = dot(sigma0q, eye(size(J, 1), size(J, 1)))
        P0 = inv(Qll)
        calib_data.statIO.sg0 = copy((dot(dot(v.T, P0), v)) / (rows - cols))
        calib_data.statIO.Exx = copy(dot(calib_data.statIO.sg0, inv(dot(dot(J.T, P0), J))))
        calib_data.statIO.varIO = copy(diag(calib_data.statIO.Exx))
        calib_data.statIO.stdIO = copy(sqrt(abs(calib_data.statIO.varIO)))
        calib_data.optimized = true
        ssc = ssout(arange(6, end()))
        calib_data.ocam_model.ss = copy(multiply(ss0, ssc.T))
        calib_data.ocam_model.xc = copy(dot(calib_data.ocam_model.xc, ssout[1]))
        calib_data.ocam_model.yc = copy(dot(calib_data.ocam_model.yc, ssout[2]))
        calib_data.ocam_model.c = ssout[3]
        calib_data.ocam_model.d = ssout[4]
        calib_data.ocam_model.e = ssout[5]
        err, stderr, MSE = reprojectpoints_quiet(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc,
                                                 calib_data.Xp_abs, calib_data.Yp_abs, M, nargout=3)
        print('Sum of squared errors:  %f\n', MSE, end='')
        MSE_old = MSE_new

    if get(gcf, 'CurrentChar') == 13:
        print('\n\nCamera model refinement interrupted', end='')
    else:
        print('\n\nCamera model optimized', end='')

    err, stderr, MSE = reprojectpoints_adv(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc,
                                           calib_data.Xp_abs, calib_data.Yp_abs, M, nargout=3)
    ss = calib_data.ocam_model.ss
    ss
    ## ================ 
    #  added code
    calib_data.errMean = err
    calib_data.errStd = stderr
    calib_data.mse = MSE
    calib_data.rms = sqrt(sum(v**2) / length(v))
    print('Root mean square[pixel]::  %f\n', calib_data.rms, end='')
    ## ================
