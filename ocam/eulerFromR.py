# this function extracts euler angles from a rotation matrix R
#     Steffen Urban email: steffen.urban@kit.edu
#     Copyright (C) 2014  Steffen Urban
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

from .libsmop import *


# this function extracts euler angles from a rotation matrix R
@function
def eulerFromR(R=None, plusPI=None):
    nargin = eulerFromR.nargin

    if nargin < 2:
        plusPI = 0

    psi2 = 0
    theta2 = 0
    phi2 = 0
    if R[3, 1] != 1 and R[3, 1] != - 1:
        theta1 = - np.arcsin(R[3, 1])
        theta2 = np.pi - theta1
        psi1 = np.arctan2(R[3, 2] / cos(theta1), R[3, 3] / cos(theta1))
        psi2 = np.arctan2(R[3, 2] / cos(theta2), R[3, 3] / cos(theta2))
        phi1 = np.arctan2(R[2, 1] / cos(theta1), R[1, 1] / cos(theta1))
        phi2 = np.arctan2(R[2, 1] / cos(theta2), R[1, 1] / cos(theta2))
        psi = psi1
        theta = theta1
        phi = phi1
    else:
        phi = 0
        if R[3, 1] == - 1:
            theta = np.pi / 2
            psi = phi + np.arctan2(R[1, 2], R[1, 3])
        else:
            theta = - np.pi / 2
            psi = - phi + np.arctan2(- R[1, 2], - R[1, 3])

    if plusPI:
        psithetaphi = concat([psi2, theta2, phi2])
    else:
        psithetaphi = concat([psi, theta, phi])

    return psithetaphi
