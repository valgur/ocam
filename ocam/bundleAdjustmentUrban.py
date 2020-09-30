from numpy import diag
from numpy.linalg import pinv

#     Steffen Urban email: steffen.urban@kit.edu
#     Copyright (C) 2014  Steffen Urban
# 
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
from .bundleErrUrban import bundleErrUrban
from .findinvpoly import findinvpoly
from .libsmop import *
from .rodrigues import rodrigues


@function
def bundleAdjustmentUrban(calib_data, robust):
    global weights
    if (robust):
        print('Starting robust non-linear refinement')
    else:
        print('Starting non-linear refinement')

    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print('\nNo linear estimate available. You must first calibrate your camera.\nClick on "Calibration"\n')
        return

    optionsLM = optimset('Display', 'off', 'TolX', 1e-05, 'TolFun', 0.0001, 'MaxIter', 100)
    if (logical_and(logical_and(isempty(calib_data.ocam_model.c), isempty(calib_data.ocam_model.d)),
                    isempty(calib_data.ocam_model.e))):
        calib_data.ocam_model.c = 1
        calib_data.ocam_model.d = 0
        calib_data.ocam_model.e = 0

    M = concat([calib_data.Xt, calib_data.Yt, zeros(size(calib_data.Xt))])
    ss0 = calib_data.ocam_model.ss
    x0 = concat([1, 1, 1, 0, 0, ones(1, size(calib_data.ocam_model.ss, 1))])
    offset = 6 + calib_data.taylor_order
    for i in calib_data.ima_proc.flat:
        R = calib_data.RRfin[:, :, i]
        R[:, 3] = np.cross(R[:, 1], R[:, 2])
        r = rodrigues(R)
        t = calib_data.RRfin[:, 3, i]
        x0 = concat([x0, r[1], r[2], r[3], t[1], t[2], t[3]])

    weights = ones(dot(dot(2, size(calib_data.Xt, 1)), length(calib_data.ima_proc)), 1)
    x0, __, vExtr, __, __, __, jacExtr = lsqnonlin(bundleErrUrban, x0, dot(- inf, ones(1, length(x0), 1)),
                                                   dot(inf, ones(1, length(x0))), optionsLM, calib_data, M, robust,
                                                   nargout=7)
    calib_data.weights = copy(weights)
    lauf = 0
    for i in calib_data.ima_proc.flat:
        RRfinOpt[:, :, i] = rodrigues(
            concat([x0[offset + 1 + lauf], x0[offset + 2 + lauf], x0[offset + 3 + lauf]]))
        RRfinOpt[:, 3, i] = concat([x0[offset + 4 + lauf], x0[offset + 5 + lauf], x0[offset + 6 + lauf]]).T
        lauf += 6

    ssc = concat([x0(arange(6, offset))])
    calib_data.ocam_model.ss = copy(multiply(ss0, ssc.T))
    calib_data.ocam_model.xc = copy(dot(calib_data.ocam_model.xc, x0[1]))
    calib_data.ocam_model.yc = copy(dot(calib_data.ocam_model.yc, x0[2]))
    calib_data.ocam_model.c = x0[3]
    calib_data.ocam_model.d = x0[4]
    calib_data.ocam_model.e = x0[5]
    ## calc standard deviation of EO
    sigma0q = 1**2
    v = copy(vExtr)
    # J = jacExtr[:,offset+1:end];
    jacExtr[:, 7] = copy([])
    J = copy(jacExtr)
    rows, cols = size(J, nargout=2)
    Qll = dot(sigma0q, eye(size[J, 1], size[J, 1]))
    P0 = inv(Qll)
    # a posteriori variance
    calib_data.statEO.sg0 = sqrt((dot(dot(v.T, P0), v)) / (rows - cols))
    # empirical covariance matrix
    Exx = pinv(dot(dot(J.T, P0), J))
    # ExxEO = inv(Jext'*Jext);
    # ExxIO = inv(Jito'*Jito);
    calib_data.statEO.Exx = copy(Exx(arange(offset, end()), arange(offset, end())))
    calib_data.statEO.varEO = copy(diag(calib_data.statEO.Exx))

    # standard deviation of ext ori parameters
    calib_data.statEO.stdEO = copy(dot(calib_data.statEO.sg0, sqrt(abs(diag(calib_data.statEO.Exx)))))
    calib_data.statIO.Exx = copy(Exx(arange(1, offset - 1), arange(1, offset - 1)))
    calib_data.statIO.varIO = copy(diag(calib_data.statIO.Exx))

    calib_data.statIO.stdIO = copy(dot(calib_data.statEO.sg0, sqrt(abs(diag(calib_data.statIO.Exx)))))
    ## rms

    calib_data.optimized = true
    calib_data.RRfin = copy(RRfinOpt)
    M = concat([calib_data.Xt, calib_data.Yt, ones(size(calib_data.Xt, 1), 1)])
    ss = calib_data.ocam_model.ss
    ss
    rms = sqrt(sum(v**2) / length(v))
    print('Root mean square[pixel]:  %f\n', rms, end='')
    calib_data.rms = copy(rms)
    calib_data.ocam_model.pol, calib_data.ocam_model.err, calib_data.ocam_model.N = findinvpoly(
        calib_data.ocam_model.ss, sqrt((calib_data.width / 2)**2 + (calib_data.height / 2)**2),
        nargout=3)
