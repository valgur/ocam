

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

from smop.libsmop import *

clear('all')
close_('all')
# !!! change those if you changed Step1_perform_test_calibrations.m !!!
nr_cams = 7
nr_tests = 5
## load scaramuzza data
for i in range(1, nr_tests + 1):
    path = sprintf('CalibData%i.mat', i)
    calibMeth[i] = load(path)

for cam in range(1, nr_cams + 1):
    for meth in range(1, nr_tests + 1):
        rms[cam, meth] = calibMeth[meth].calib_data[cam].rms

figure
bar(rms)
title('root mean square error')
xlabel('data set')
ylabel('root mean square error [pixel]')
legend('ocam standard', 'ocam subpixel', 'urban', 'urban subpixel', 'urban subpixel robust')
for cam in range(1, nr_cams + 1):
    for meth in range(1, nr_tests + 1):
        runtime[cam, meth] = calibMeth[meth].calib_data[cam].runtime

figure
bar(runtime)
title('runtime')
xlabel('data set')
ylabel('time [s]')
legend('ocam standard', 'ocam subpixel', 'urban', 'urban subpixel', 'urban subpixel robust')
badIdx = copy([])
for cam in range(1, nr_cams + 1):
    for meth in range(1, 2 + 1):
        stdEOangle = copy([])
        stdEOpos = copy([])
        for imgs in calibMeth[meth].calib_data[cam].ima_proc.flat:
            try:
                if (logical_not(isempty(calibMeth[meth].calib_data[cam].statEO[imgs].stdEO))):
                    stdEOangle = concat([stdEOangle, dot(1000.0, concat(
                        [calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[1],
                         calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[2],
                         calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[3]]))])
                    stdEOpos = concat([stdEOpos, concat([calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[4],
                                                         calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[5],
                                                         calibMeth[meth].calib_data[cam].statEO[imgs].stdEO[6]])])
                    badIdx = concat([badIdx, imgs])
            finally:
                pass
        stdEOangleA[cam, meth] = mean(stdEOangle)
        stdEOposA[cam, meth] = mean(stdEOpos)
        badIdx = copy([])

lauf = 0
for cam in range(1, nr_cams + 1):
    for meth in arange(3, nr_tests).flat:
        stdEOangle = copy([])
        stdEOpos = copy([])
        for imgs in arange(1, length(calibMeth[meth].calib_data[cam].ima_proc)).flat:
            stdEOangle = concat([stdEOangle, dot(1000.0, concat([calibMeth[meth].calib_data[cam].statEO.stdEO(1 + lauf),
                                                                 calibMeth[meth].calib_data[cam].statEO.stdEO(2 + lauf),
                                                                 calibMeth[meth].calib_data[cam].statEO.stdEO(
                                                                     3 + lauf)]))])
            stdEOpos = concat([stdEOpos, concat([calibMeth[meth].calib_data[cam].statEO.stdEO(4 + lauf),
                                                 calibMeth[meth].calib_data[cam].statEO.stdEO(5 + lauf),
                                                 calibMeth[meth].calib_data[cam].statEO.stdEO(6 + lauf)])])
            lauf += 6
        lauf = 0
        stdEOangleA[cam, meth] = median(stdEOangle)
        stdEOposA[cam, meth] = median(stdEOpos)

figure
bar(stdEOangleA)
title('mean standard deviation of camera orientations')
xlabel('data set')
ylabel('mean standard deviation[mrad]')
legend('ocam standard', 'ocam subpixel', 'urban', 'urban subpixel', 'urban subpixel robust')
figure
bar(stdEOposA)
title('mean standard deviation of camera positions')
xlabel('data set')
ylabel('standard deviation [mm]')
legend('ocam standard', 'ocam subpixel', 'urban', 'urban subpixel', 'urban subpixel robust')
## IO, c and a0
lauf = 0
for cam in range(1, nr_cams + 1):
    for meth in range(1, nr_tests + 1):
        std_cde[cam, meth] = calibMeth[meth].calib_data[cam].statIO.stdIO[3]
        std_a0[cam, meth] = calibMeth[meth].calib_data[cam].statIO.stdIO[6]
