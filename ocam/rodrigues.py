# RODRIGUES	Transform rotation matrix into rotation vector and viceversa.
#
#		Syntax:  [OUT]=RODRIGUES(IN)
# 		If IN is a 3x3 rotation matrix then OUT is the
#		corresponding 3x1 rotation vector
# 		if IN is a rotation 3-vector then OUT is the 
#		corresponding 3x3 rotation matrix
#
##
##		Copyright (c) March 1993 -- Pietro Perona
##		California Institute of Technology
##
#
## ALL CHECKED BY JEAN-YVES BOUGUET, October 1995.
## FOR ALL JACOBIAN MATRICES !!! LOOK AT THE TEST AT THE END !!
#
## BUG when norm(om)=np.pi fixed -- April 6th, 1997;
## Jean-Yves Bouguet
#
## Add projection of the 3x3 matrix onto the set of special ortogonal matrices SO[3] by SVD -- February 7th, 2003;
## Jean-Yves Bouguet

from .libsmop import *


@function
def rodrigues(in_):
    m, n = size(in_, nargout=2)
    # bigeps = 10e+4*eps;
    bigeps = dot(1e+21, eps)
    if logical_or((logical_and((m == 1), (n == 3))), (logical_and((m == 3), (n == 1)))):
        theta = norm(in_)
        if theta < eps:
            R = eye[3]
            dRdin = concat(
                [[0, 0, 0], [0, 0, 1], [0, - 1, 0], [0, 0, - 1], [0, 0, 0], [1, 0, 0], [0, 1, 0], [- 1, 0, 0],
                 [0, 0, 0]])
        else:
            if n == length(in_):
                in_ = in_.T
            # m3 = [in,theta]
            dm3din = concat([[eye[3]], [in_.T / theta]])
            omega = in_ / theta
            dm2dm3 = concat([[eye[3] / theta, - in_ / theta**2], [zeros(1, 3), 1]])
            alpha = np.cos(theta)
            beta = sin(theta)
            gamma = 1 - np.cos(theta)
            omegav = concat([[concat([0, - omega[3], omega[2]])], [concat([omega[3], 0, - omega[1]])],
                             [concat([- omega[2], omega[1], 0])]])
            A = dot(omega, omega.T)
            dm1dm2 = zeros(21, 4)
            dm1dm2[1, 4] = - sin(theta)
            dm1dm2[2, 4] = np.cos(theta)
            dm1dm2[3, 4] = sin(theta)
            dm1dm2[arange(4, 12), arange(1, 3)] = concat(
                [[0, 0, 0, 0, 0, 1, 0, - 1, 0], [0, 0, - 1, 0, 0, 0, 1, 0, 0], [0, 1, 0, - 1, 0, 0, 0, 0, 0]]).T
            w1 = omega[1]
            w2 = omega[2]
            w3 = omega[3]
            dm1dm2[arange(13, 21), 1] = concat([[dot(2, w1)], [w2], [w3], [w2], [0], [0], [w3], [0], [0]])
            dm1dm2[arange(13, 21), 2] = concat([[0], [w1], [0], [w1], [dot(2, w2)], [w3], [0], [w3], [0]])
            dm1dm2[arange(13, 21), 3] = concat([[0], [0], [w1], [0], [0], [w2], [w1], [w2], [dot(2, w3)]])
            R = dot(eye[3], alpha) + dot(omegav, beta) + dot(A, gamma)
            dRdm1 = zeros(9, 21)
            dRdm1[concat([1, 5, 9]), 1] = ones(3, 1)
            dRdm1[:, 2] = ravel(omegav)
            dRdm1[:, arange(4, 12)] = dot(beta, eye[9])
            dRdm1[:, 3] = ravel(A)
            dRdm1[:, arange(13, 21)] = dot(gamma, eye[9])
            dRdin = dot(dot(dot(dRdm1, dm1dm2), dm2dm3), dm3din)
        out = copy(R)
        dout = copy(dRdin)
    else:
        if (logical_and(logical_and(logical_and((m == n), (m == 3)), (norm(dot(in_.T, in_) - eye[3]) < bigeps)),
                        (abs(det(in_) - 1) < bigeps))):
            R = copy(in_)
            U, S, V = svd(R, nargout=3)
            R = dot(U, V.T)
            tr = (trace(R) - 1) / 2
            dtrdR = concat([1, 0, 0, 0, 1, 0, 0, 0, 1]) / 2
            theta = real(acos(tr))
            if sin(theta) >= 1e-05:
                dthetadtr = - 1 / sqrt(1 - tr**2)
                dthetadR = dot(dthetadtr, dtrdR)
                vth = 1 / (dot(2, sin(theta)))
                dvthdtheta = dot(- vth, np.cos(theta)) / sin(theta)
                dvar1dtheta = concat([[dvthdtheta], [1]])
                dvar1dR = dot(dvar1dtheta, dthetadR)
                om1 = concat([R(3, 2) - R(2, 3), R(1, 3) - R(3, 1), R(2, 1) - R(1, 2)]).T
                dom1dR = concat(
                    [[0, 0, 0, 0, 0, 1, 0, - 1, 0], [0, 0, - 1, 0, 0, 0, 1, 0, 0], [0, 1, 0, - 1, 0, 0, 0, 0, 0]])
                dvardR = concat([[dom1dR], [dvar1dR]])
                om = dot(vth, om1)
                domdvar = concat([dot(vth, eye[3]), om1, zeros(3, 1)])
                dthetadvar = concat([0, 0, 0, 0, 1])
                dvar2dvar = concat([[domdvar], [dthetadvar]])
                out = dot(om, theta)
                domegadvar2 = concat([dot(theta, eye[3]), om])
                dout = dot(dot(domegadvar2, dvar2dvar), dvardR)
            else:
                if tr > 0:
                    out = concat([0, 0, 0]).T
                    dout = concat([[0, 0, 0, 0, 0, 1 / 2, 0, - 1 / 2, 0], [0, 0, - 1 / 2, 0, 0, 0, 1 / 2, 0, 0],
                                   [0, 1 / 2, 0, - 1 / 2, 0, 0, 0, 0, 0]])
                else:
                    out = dot(theta, (
                        multiply(sqrt((diag(R) + 1) / 2), concat([[1], [dot(2, (R(1, arange(2, 3)) >= 0).T) - 1]]))))
                    if nargout > 1:
                        print('WARNING!!!! Jacobian domdR undefined!!!')
                        dout = dot(np.nan, ones(3, 9))
        else:
            error('Neither a rotation matrix nor a rotation vector were provided')

    return out, dout


def test_rodrigues():
    ## test of the Jacobians:

    #### TEST OF dRdom:
    om = randn(3, 1)
    dom = randn(3, 1) / 1000000
    R1, dR1 = rodrigues(om, nargout=2)
    R2 = rodrigues(om + dom)
    R2a = R1 + reshape(dot(dR1, dom), 3, 3)
    gain = norm(R2 - R1) / norm(R2 - R2a)
    ### TEST OF dOmdR:
    om = randn(3, 1)
    R = rodrigues(om)
    dom = randn(3, 1) / 10000
    dR = rodrigues(om + dom) - R
    omc, domdR = rodrigues(R, nargout=2)
    om2 = rodrigues(R + dR)
    om_app = omc + dot(domdR, ravel(dR))
    gain = norm(om2 - omc) / norm(om2 - om_app)
    ### OTHER BUG: (FIXED NOW!!!)

    omu = randn(3, 1)
    omu /= norm(omu)
    om = dot(np.pi, omu)
    R, dR = rodrigues(om, nargout=2)
    om2 = rodrigues(R)
    concat([om, om2])
    ### NORMAL OPERATION

    om = randn(3, 1)
    R, dR = rodrigues(om, nargout=2)
    om2 = rodrigues(R)
    concat([om, om2])
