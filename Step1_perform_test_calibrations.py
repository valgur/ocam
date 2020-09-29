# Steffen Urban email: steffen.urban@kit.edu
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

from smop.libsmop import *

from ocam.bundleAdjustmentUrban import bundleAdjustmentUrban
from ocam.calib_data import CalibData
from ocam.calibration import calibration
from ocam.check_active_images import check_active_images
from ocam.check_directory import check_directory
from ocam.findcenter import findcenter
from ocam.get_checkerboard_corners import get_checkerboard_corners
from ocam.optimizefunction import optimizefunction

# test cases: [use_urban, use_subpixel, robust]
test_cases = cellarray([
    concat([false, false, false]),
    concat([false, true, false]),
    concat([true, false, false]),
    concat([true, true, false]),
    concat([true, true, true])
])
nr_tests = size(test_cases, 2)
# filenames
base_names = cellarray(['Fisheye1_', 'Fisheye2_', 'GOPR', 'MiniOmni', 'VMRImage', 'Ladybug', 'KaidanOmni'])

# nr of squares
squares = cellarray([
    concat([5, 7]),
    concat([5, 7]),
    concat([5, 7]),
    concat([6, 10]),
    concat([5, 6]),
    concat([4, 7]),
    concat([6, 10])
])
# size of squares in [mm]
sizes = cellarray([32.5, 117, 117, 30, 30, 30, 30])
# change polynomial degree if necessary
polDegree = cellarray([4, 4, 4, 4, 4, 4, 4])
corners_already_extracted = 1

## this loop automatically calls all relevant calibration functions
# loop over test scenarios
for b in arange(1, size(test_cases, 2)).flat:
    # loop over data sets
    calib_data = copy([])
    for idx in arange(1, size(base_names, 2)).flat:
        # set bool variables
        use_urban = test_cases[b](1)
        use_subpixel = test_cases[b](2)
        calib_data[idx] = CalibData()
        calib_data[idx].calib_name = base_names[idx]
        Nima_valid = check_directory(calib_data[idx])
        ima_read_calib(calib_data[idx])
        calib_data[idx].calibrated = 0
        check_active_images(calib_data[idx])
        calib_data[idx].taylor_order = polDegree[idx]
        calib_data[idx].taylor_order_default = polDegree[idx]
        use_video_mode = 1
        use_corner_find = 0
        calib_data[idx].Xt = copy([])
        calib_data[idx].Yt = copy([])
        calib_data[idx].ocam_model.xc = round(calib_data[idx].ocam_model.height / 2)
        calib_data[idx].ocam_model.yc = round(calib_data[idx].ocam_model.width / 2)
        calib_data[idx].dX = sizes[idx]
        calib_data[idx].dY = sizes[idx]
        calib_data[idx].n_sq_x = squares[idx](1)
        calib_data[idx].n_sq_y = squares[idx](2)
        for i in arange(0, calib_data[idx].n_sq_x).flat:
            for j in arange(0, calib_data[idx].n_sq_y).flat:
                calib_data[idx].Yt = concat([[calib_data[idx].Yt], [dot(j, calib_data[idx].dY)]])
                calib_data[idx].Xt = concat([[calib_data[idx].Xt], [dot(i, calib_data[idx].dX)]])
        ima_numbers = arange(1, calib_data[idx].n_imgs)
        count = 0
        for kk in ima_numbers.flat:
            callBack, x, y = get_checkerboard_corners(kk, use_subpixel, calib_data[idx], nargout=3)
            if callBack == 1:
                count += 1
                calib_data[idx].Xp_abs[:, :, kk] = x
                calib_data[idx].Yp_abs[:, :, kk] = y
                calib_data[idx].active_images[kk] = 1
                calib_data[idx].ima_proc = copy(sort(concat([calib_data[idx].ima_proc, kk])))
            calib_data[idx].I[kk] = 1
        # perform linear calibration
        # this step is equal to all methods
        calibration(calib_data[idx])
        if use_urban:
            # bundle adjustment, can be robust
            bundleAdjustmentUrban(calib_data[idx], test_cases[b](3))
        else:
            # ocam_calib method
            findcenter(calib_data[idx])
            optimizefunction(calib_data[idx])
        calib_data[idx].runtime = toc
    # save results
    path = f'CalibData{b:d}.mat'
    save(path, 'calib_data')
