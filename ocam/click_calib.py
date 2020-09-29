from . import click_ima_calib, click_ima_calib_rufli
from .calib_data import CalibData
from .get_checkerboard_corners import get_checkerboard_corners
from .libsmop import *


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
# 04.03.2014 by Steffen Urban
# this is a modified file from
# Davide Scaramuzzas Toolbox OcamCalib
# original filename: click_calib.m

@function
def click_calib(calib_data: CalibData, ima_numbers=None, use_video_mode=True, use_automatic=True):
    if ima_numbers is None:
        ima_numbers = arange(1, len(calib_data.ima_proc))

    # Arranging the pixel of the world
    calib_data.Xt = copy([])
    calib_data.Yt = copy([])
    for i in arange(0, calib_data.n_sq_x).flat:
        for j in arange(0, calib_data.n_sq_y).flat:
            calib_data.Yt = concat([
                [calib_data.Yt],
                [j * calib_data.dY]
            ])
            calib_data.Xt = concat([
                [calib_data.Xt],
                [i * calib_data.dX]
            ])

    use_corner_find = True
    if use_video_mode:
        count = 0
        for kk in ima_numbers:
            callBack, x, y = get_checkerboard_corners(kk, use_corner_find, calib_data, nargout=3)
            if callBack == 1:
                count += 1
                calib_data.Xp_abs[:, :, kk] = x
                calib_data.Yp_abs[:, :, kk] = y
                calib_data.active_images[kk] = 1
                calib_data.ima_proc = matlabarray(sort(concat([calib_data.ima_proc, kk])))
    else:
        for kk in ima_numbers:
            if not use_automatic:
                click_ima_calib(kk, use_corner_find, calib_data)
            else:
                click_ima_calib_rufli(kk, use_corner_find, calib_data)
            calib_data.active_images[kk] = True
            calib_data.ima_proc = matlabarray(sort(concat([calib_data.ima_proc, kk])))
