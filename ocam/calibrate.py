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

from numpy import sign
from numpy.linalg import svd, norm, pinv

from .libsmop import *


@function
def calibrate(Xt=None, Yt=None, Xp_abs=None, Yp_abs=None, xc=None, yc=None, taylor_order_default=None, ima_proc=None,
              *args, **kwargs):
    Yp = Yp_abs - yc
    Xp = Xp_abs - xc
    for kk in ima_proc.flat:
        Ypt = Yp[:, :, kk]
        Xpt = Xp[:, :, kk]
        A = concat([multiply(Xt, Ypt), multiply(Yt, Ypt), multiply(- Xt, Xpt), multiply(- Yt, Xpt), Ypt, - Xpt])
        U, S, V = svd(A, nargout=3)
        R11 = V(1, end())
        R12 = V(2, end())
        R21 = V(3, end())
        R22 = V(4, end())
        T1 = V(5, end())
        T2 = V(6, end())
        AA = ((dot(R11, R12)) + (dot(R21, R22)))**2
        BB = R11**2 + R21**2
        CC = R12**2 + R22**2
        R32_2 = roots(concat([1, CC - BB, - AA]))
        R32_2 = R32_2(find(R32_2 >= 0))
        R31 = copy([])
        R32 = copy([])
        sg = concat([1, - 1])
        for i in arange(1, length(R32_2)).flat:
            for j in arange(1, 2).flat:
                sqrtR32_2 = dot(sg(j), sqrt(R32_2(i)))
                R32 = concat([[R32], [sqrtR32_2]])
                if R32_2 == 0:
                    R31 = concat([[R31], [sqrt(CC - BB)]])
                    R31 = concat([[R31], [-sqrt(CC - BB)]])
                    R32 = concat([[R32], [sqrtR32_2]])
                else:
                    R31 = concat([[R31], [numpy.linalg.solve(- sqrtR32_2, (dot(R11, R12) + dot(R21, R22)))]])
        RR = zeros(3, 3, dot(length(R32), 2))
        count = 0
        for i1 in arange(1, length(R32)).flat:
            for i2 in arange(1, 2).flat:
                count += 1
                Lb = numpy.linalg.solve(sqrt(R11**2 + R21**2 + R31(i1)**2), 1)
                RR[:, :, count] = dot(dot(sg(i2), Lb),
                                      concat([[R11, R12, T1], [R21, R22, T2], [R31(i1), R32(i1), 0]]))
        RR1 = copy([])
        minRR = np.inf
        minRR_ind = - 1
        for min_count in arange(1, size(RR, 3)).flat:
            if (norm(concat([[RR(1, 3, min_count)], [RR(2, 3, min_count)]]) - concat([[Xpt[1]], [Ypt[1]]])) < minRR):
                minRR = norm(concat([[RR(1, 3, min_count)], [RR(2, 3, min_count)]]) - concat([[Xpt[1]], [Ypt[1]]]))
                minRR_ind = min_count
        if minRR_ind != - 1:
            count2 = 0
            for count in arange(1, size(RR, 3)).flat:
                if (sign(RR(1, 3, count)) == logical_and(sign(RR(1, 3, minRR_ind)), sign(RR(2, 3, count))) == sign(
                        RR(2, 3, minRR_ind))):
                    count2 += 1
                    RR1[:, :, count2] = RR[:, :, count]
        if isempty(RR1):
            RRfin = 0
            ss = 0
            return RRfin, ss
        #    figure[1]; imagesc(aa(:,:,:,counter)); set(h,'name',filename);
        nm = plot_RR(RR1, Xt, Yt, Xpt, Ypt, 0)
        RRdef = RR1[:, :, nm]
        RRfin[:, :, kk] = RRdef

    # Run program finding mirror and translation parameters
    RRfin, ss = omni_find_parameters_fun(Xt, Yt, Xp_abs, Yp_abs, xc, yc, RRfin, taylor_order_default, ima_proc,
                                         nargout=2)
    return RRfin, ss


@function
def omni_find_parameters_fun(Xt=None, Yt=None, Xp_abs=None, Yp_abs=None, xc=None, yc=None, RRfin=None,
                             taylor_order=None, ima_proc=None):
    Yp = Yp_abs - yc
    Xp = Xp_abs - xc
    # obrand_start
    min_order = 4
    range_ = 1000
    if taylor_order <= min_order:
        PP = copy([])
        QQ = copy([])
        mono = copy([])
        count = 0
        for i in ima_proc.flat:
            count += 1
            RRdef = RRfin[:, :, i]
            R11 = RRdef(1, 1)
            R21 = RRdef(2, 1)
            R31 = RRdef(3, 1)
            R12 = RRdef(1, 2)
            R22 = RRdef(2, 2)
            R32 = RRdef(3, 2)
            T1 = RRdef(1, 3)
            T2 = RRdef(2, 3)
            Xpt = Xp[:, :, i]
            Ypt = Yp[:, :, i]
            MA = multiply(R21, Xt) + multiply(R22, Yt) + T2
            MB = multiply(Ypt, (multiply(R31, Xt) + multiply(R32, Yt)))
            MC = multiply(R11, Xt) + multiply(R12, Yt) + T1
            MD = multiply(Xpt, (multiply(R31, Xt) + multiply(R32, Yt)))
            rho = copy([])
            for j in arange(2, taylor_order).flat:
                rho[:, :, j] = (sqrt(Xpt**2 + Ypt**2))**j
            PP1 = concat([[MA], [MC]])
            for j in arange(2, taylor_order).flat:
                PP1 = concat([PP1, concat(
                    [[multiply(MA, rho[:, :, j])], [multiply(MC, rho[:, :, j])]])])
            PP = concat(
                [[PP, zeros(size(PP, 1), 1)], [PP1, zeros(size(PP1, 1), count - 1), concat([[- Ypt], [- Xpt]])]])
            QQ = concat([[QQ], [MB], [MD]])
        if (license('checkout', 'Optimization_Toolbox') != 1):
            s = dot(pinv(PP), QQ)
            ss = s(arange(1, taylor_order))
            count = 0
            for j in ima_proc.flat:
                count += 1
                RRfin[3, 3, j] = s(length(ss) + count)
            ss = concat([[ss[1]], [0], [ss(arange(2, end()))]])
        else:
            diff_order = 2
            mono_rho = copy([])
            for j in arange(1, taylor_order).flat:
                mono_rho[:, j] = concat([arange(1, max(squeeze(max(Xp_abs))))])**j
            for diff_k in arange(1, min(diff_order, taylor_order)).flat:
                mono1 = copy([])
                for j in arange(1, diff_k).flat:
                    mono1 = concat([mono1, zeros(size(mono_rho, 1), 1)])
                mono1 = concat([mono1, dot(factorial(diff_k), ones(size(mono_rho, 1), 1))])
                for j in arange(diff_k + 1, taylor_order).flat:
                    mono1 = concat([mono1, dot(factorial(j) / factorial(j - diff_k), mono_rho[:, j - diff_k])])
                mono = concat([[mono], [mono1[:, 1], mono1[:, arange(3, end())]]])
            if isempty(mono):
                mono = copy([])
                mono_const = copy([])
            else:
                mono = - concat([mono, zeros(size(mono, 1), size(PP, 2) - size(mono, 2))])
                mono_const = zeros(size(mono, 1), 1)
            options = optimset('Display', 'off')
            warning('off', 'all')
            s_o, ob_resnorm, ob_residual, ob_exitflag, ob_output, ob_lambda = lsqlin(PP, QQ, mono, mono_const, [], [],
                                                                                     [], [], [], options, nargout=6)
            warning('on', 'all')
            ss_o = s_o(arange(1, taylor_order))
            count = 0
            for j in ima_proc.flat:
                count += 1
                RRfin[3, 3, j] = s_o(length(ss_o) + count)
            ss = concat([[ss_o[1]], [0], [ss_o(arange(2, end()))]])
    else:
        # obrand calculate higher order polynomials iteratively
        x0 = zeros(1, length(ima_proc) + 2)
        lb = dot(- np.inf, ones(size(x0)))
        ub = dot(np.inf, ones(size(x0)))
        max_order = taylor_order
        for taylor_order in arange(min_order, max_order).flat:
            PP = copy([])
            QQ = copy([])
            count = 0
            for i in ima_proc.flat:
                count += 1
                RRdef = RRfin[:, :, i]
                R11 = RRdef(1, 1)
                R21 = RRdef(2, 1)
                R31 = RRdef(3, 1)
                R12 = RRdef(1, 2)
                R22 = RRdef(2, 2)
                R32 = RRdef(3, 2)
                T1 = RRdef(1, 3)
                T2 = RRdef(2, 3)
                Xpt = Xp[:, :, i]
                Ypt = Yp[:, :, i]
                MA = multiply(R21, Xt) + multiply(R22, Yt) + T2
                MB = multiply(Ypt, (multiply(R31, Xt) + multiply(R32, Yt)))
                MC = multiply(R11, Xt) + multiply(R12, Yt) + T1
                MD = multiply(Xpt, (multiply(R31, Xt) + multiply(R32, Yt)))
                rho = copy([])
                for j in arange(2, taylor_order).flat:
                    rho[:, :, j] = (sqrt(Xpt**2 + Ypt**2))**j
                PP1 = concat([[MA], [MC]])
                for j in arange(2, taylor_order).flat:
                    PP1 = concat([PP1, concat(
                        [[multiply(MA, rho[:, :, j])], [multiply(MC, rho[:, :, j])]])])
                PP = concat(
                    [[PP, zeros(size(PP, 1), 1)], [PP1, zeros(size(PP1, 1), count - 1), concat([[- Ypt], [- Xpt]])]])
                QQ = concat([[QQ], [MB], [MD]])
            if (license('checkout', 'Optimization_Toolbox') != 1):
                s = dot(pinv(PP), QQ)
                ss = s(arange(1, taylor_order))
            else:
                options = optimset('Display', 'off')
                warning('off', 'all')
                s, ob_resnorm, ob_residual, ob_exitflag, ob_output, ob_lambda = lsqlin(PP, QQ, [], [], [], [], lb, ub,
                                                                                       x0, options, nargout=6)
                warning('on', 'all')
                ss = s(arange(1, taylor_order))
            x0 = concat([[s(arange(1, taylor_order))], [0], [s(arange(taylor_order + 1, end()))]])
            lb = concat([[s(arange(1, taylor_order)) - abs(dot(range_, s(arange(1, taylor_order))))], [- np.inf],
                         [dot(- np.inf, ones(size(s(arange(taylor_order + 1, end())))))]])
            ub = concat([[s(arange(1, taylor_order)) + abs(dot(range_, s(arange(1, taylor_order))))], [np.inf],
                         [dot(np.inf, ones(size(s(arange(taylor_order + 1, end())))))]])
            eq_bounds = find(lb == ub)
            lb[eq_bounds] = lb(eq_bounds) - eps(lb(eq_bounds))
            ub[eq_bounds] = ub(eq_bounds) + eps(ub(eq_bounds))
        count = 0
        for j in ima_proc.flat:
            count += 1
            RRfin[3, 3, j] = s(length(ss) + count)
        ss = concat([[ss[1]], [0], [ss(arange(2, end()))]])

    # obrand_end
    return RRfin, ss


@function
def plot_RR(RR=None, Xt=None, Yt=None, Xpt=None, Ypt=None, figure_number=None):
    if figure_number > 0:
        figure(figure_number)

    for i in arange(1, size(RR, 3)).flat:
        RRdef = RR[:, :, i]
        R11 = RRdef(1, 1)
        R21 = RRdef(2, 1)
        R31 = RRdef(3, 1)
        R12 = RRdef(1, 2)
        R22 = RRdef(2, 2)
        R32 = RRdef(3, 2)
        T1 = RRdef(1, 3)
        T2 = RRdef(2, 3)
        MA = multiply(R21, Xt) + multiply(R22, Yt) + T2
        MB = multiply(Ypt, (multiply(R31, Xt) + multiply(R32, Yt)))
        MC = multiply(R11, Xt) + multiply(R12, Yt) + T1
        MD = multiply(Xpt, (multiply(R31, Xt) + multiply(R32, Yt)))
        rho = sqrt(Xpt**2 + Ypt**2)
        rho2 = (Xpt**2 + Ypt**2)
        PP1 = concat([[MA, multiply(MA, rho), multiply(MA, rho2)], [MC, multiply(MC, rho), multiply(MC, rho2)]])
        PP = concat([PP1, concat([[- Ypt], [- Xpt]])])
        QQ = concat([[MB], [MD]])
        s = dot(pinv(PP), QQ)
        ss = s(arange(1, 3))
        if figure_number > 0:
            subplot(1, size(RR, 3), i)
            plot(arange(0, 620), polyval(concat([ss[3], ss[2], ss[1]]), concat([arange(0, 620)])))
            grid
            axis('equal')
        if ss(end()) >= 0:
            index = i
    return index
