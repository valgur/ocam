#     Steffen Urban email: steffen.urban@kit.edu
#     Copyright (C) 2014 Steffen Urban, 2009 Davide Scaramuzza
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
#
# 04.03.2014 by Steffen Urban
# this is a modified file from
# Davide Scaramuzzas Toolbox OcamCalib
# filename: get_checkerboard_corners.m
# Code was changed to force subpixel corners

import os
import shutil
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np

from . import libsmop
from .calib_data import CalibData
from .cornerfinder import cornerfinder
from .draw_axes import draw_axes


def get_checkerboard_corners(calib_data: CalibData, kk: int, refine_corners: bool = True, visualize: bool = False):
    n_sq_x, n_sq_y = calib_data.n_sq_x, calib_data.n_sq_y
    img = calib_data.read_image(kk) if refine_corners or visualize else None

    print(f'Processing image {calib_data.imgs[kk - 1]}...')

    corner_info, cornersX, cornersY = call_find_corners(calib_data.imgs[kk - 1], n_sq_x, n_sq_y)
    if corner_info is None:
        return None, None

    # cornersX, cornersY = _fix_ambiguous_corners(corner_info, cornersX, cornersY)
    if refine_corners:
        cornersX, cornersY = _refine_corners(cornersX, cornersY, img)
    # _fix_orientation_ambiguity()
    # _fix_orientation_ambiguity2()
    x, y = _flatten_corner_matrices(cornersX, cornersY)
    x, y = _fix_numbering_direction(x, y, n_sq_x, n_sq_y)

    if n_sq_x != n_sq_y:
        x, y = _fix_nonsquare_board_ambiguity(x, y, n_sq_x, n_sq_y)

    if visualize:
        visualize_corners(x, y, img, f'Image {kk}', n_sq_y)

    print('Done')

    # convert to Matlab format for compatibility
    x[x >= 0] += 1
    y[x >= 0] += 1
    x = libsmop.copy(x)
    y = libsmop.copy(y)

    # return swapped x and y
    return y.T, x.T


def call_find_corners(img_file: str, n_sq_x: int, n_sq_y: int) -> (np.ndarray, np.ndarray, np.ndarray):
    img_file = os.path.abspath(img_file)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    executable = os.path.join(script_dir, 'autoCornerFinder/FindCorners')
    cur_dir = os.curdir
    work_dir = tempfile.mkdtemp()
    try:
        os.chdir(work_dir)
        os.mkdir('cToMatlab')

        # Tell the automatic corner extractor, which image file to process
        with open('pictures.txt', 'w') as f:
            f.write(img_file)

        # Call the automatic corner extraction algorithm
        # Width and height in the algorithm are defined differently than
        # in this toolbox: Its the number of inernal corners instead of
        # internal quadrangles.
        # -m specifies the number of corners the algorithm has to find, before it
        # terminates
        # If retcode = -1   -> An error occured, automatic corner finding is aborted
        #            =  0   -> Not enough corners have been found, add some manually
        #            =  1   -> Enough corners have been found for calibration

        retcode = subprocess.call([
            executable,
            '-w', str(n_sq_x + 1),
            '-h', str(n_sq_y + 1),
            '-m', str((n_sq_x + 1) * (n_sq_y + 1)),
            'pictures.txt'
        ])

        # Do error checking
        if retcode == 0:
            # Image omitted -- Not all corners found
            return None, None, None
        elif retcode != 1:
            with open('cToMatlab/error.txt') as f:
                line = f.read()
            raise RuntimeError('Image omitted -- During corner finding an error occurred: ' + line)

        # Open the files with the found corners
        cornerInfo = np.loadtxt('cToMatlab/cornerInfo.txt', dtype=int)
        cornersX = np.loadtxt('cToMatlab/cornersX.txt')
        cornersY = np.loadtxt('cToMatlab/cornersY.txt')
    finally:
        os.chdir(cur_dir)
        shutil.rmtree(work_dir)
    return cornerInfo, cornersX, cornersY


def _fix_ambiguous_corners(cornerInfo, cornersX, cornersY):
    # If there are entire rows of zeros at the end of the file, this means that
    # the cornerfinder is not entirely sure about the location of the board.
    # Then eventually additional corners are needed in the front.
    if np.any(~np.any(cornersX >= 0, axis=1)):
        raise NotImplementedError
        # Add at least one row! This is because the columns are maybe ambiguous, 
        # and then we don't know whether they belong to the current row or to the
        # next one...
        # if cornerInfo[0] - iSave >= 0:
        #     cornersX = concat([
        #         -1 * ones(cornerInfo[0] - iSave + 1, cornerInfo[1]),
        #         cornersX
        #     ])
        #     cornersY = concat([
        #         -1 * ones(cornerInfo[0] - iSave + 1, cornerInfo[1]),
        #         cornersY
        #     ])
    return cornersX, cornersY


def _fix_orientation_ambiguity():
    pass
    # #ORIENTATION AMBIGUITY?
    # ###########################################################################
    # #Decide whether the position and orientation of the checkerboard
    # #can be unambiguously determined
    # if min((n_sq_y+1), (n_sq_x+1)) ~= min(cornerInfo[0], cornerInfo[1])
    #     if flagStart == true
    #         #DISPLAY AN ERROR AND RETURN
    #     end
    #     #There is some ambiguity, get rid of it
    #     figure[2];
    #     imagesc(I);
    #     colormap(gray);
    #     set(2,'color',[1 1 1]);
    #     title({['Image ' str(kk)]});
    #     h = get(gca, 'title');
    #     set(h, fontweight='bold')
    #     axis([min_x max_x min_y max_y]);
    #     
    # 
    #     figure[2]; hold on;
    #     i = startingCorner[1]
    #     j = startingCorner[2]
    #     #Plot the starting corner and its neighbor
    #     plot( cornersX(i,j),cornersY(i,j),'+','color','red','linewidth',2);
    #     plot( cornersX(i,j+1),cornersY(i,j+1),'+','color','red','linewidth',2);
    #     text( cornersX(i,j)+3,cornersY(i,j)+3,str[0] )
    #     text( cornersX(i,j+1)+3,cornersY(i,j+1)+3,str[1] )
    #     set(findobj('type', 'text'), color='red'); 
    #     hold off;
    #       
    #     #Get user input on behalf of the board orientation
    #     numCornersDirZeroOne = input(['The automatic corner finder was not able to decide how the pattern is oriented. Please indicate the number of corners in direction of corners [0 -> 1]: ']);
    #     #We look in row direction
    #     deltaCols = cornerInfo[1] - numCornersDirZeroOne;                 
    # end


def _refine_corners(cornersX: np.ndarray, cornersY: np.ndarray, img: np.array) -> (np.ndarray, np.ndarray):
    # Apply the corner finder
    # added by steffen urban
    rows, cols = cornersX.shape
    height, width = img.shape
    wintx = int(np.ceil(height / 100))
    winty = int(np.ceil(height / 100))
    for i in range(rows):
        for j in range(cols):
            if cornersX[i, j] >= 0:
                x, y = cornerfinder(np.array([[cornersX[i, j], cornersY[i, j]]]).T, img, winty, wintx)[0]
                cornersX[i, j] = x
                cornersY[i, j] = y
    return cornersX, cornersY


def _flatten_corner_matrices(cornersX: np.ndarray, cornersY: np.ndarray) -> (np.ndarray, np.ndarray):
    # Save all corners in two arrays for further processing by other functions
    rows, cols = cornersX.shape
    num_corners = rows * cols
    # Find the first existing corner
    min_i = None
    min_j = None
    for i in range(rows):
        for j in range(cols):
            if cornersX[i, j] != -1 and min_i is None:
                min_i = i
                min_j = j
                break
    # Start with the smallest found corner and then append the larger ones
    x = []
    y = []
    iteration = 0
    while True:
        # Iterators
        i = (min_i + 1) + iteration // cols - 1
        j = ((min_j + iteration) % cols)
        iteration += 1
        x.append(float(cornersX[i, j]))
        y.append(float(cornersY[i, j]))
        if iteration >= num_corners or i >= rows:
            break
    return np.array(x), np.array(y)


def _fix_orientation_ambiguity2():
    pass
    #     #DOES NOT WORK RELIABLY!
    #     #3. In case of a n x m board, where n is even and m is odd (or vice 
    #     #versa) there still exists a 180 degree ambiguity. This could be get 
    #     #rid off here.
    #     #We define the starting corner as the corner where "a black checker is
    #     #outernmost". Only proceed if one dimension is odd and the other even!
    #     if (mod(n_cor_min,2) ~= 0 & mod(n_cor_max,2) == 0) | (mod(n_cor_min,2) == 0 & mod(n_cor_max,2) ~= 0)
    #         #Determine the intensity at the given locations in the image
    #         [dummy,lengthl] = size(x);
    #         intens1 = I( floor((y(1,1)+y(1,n_cor_max+2))/2), floor((x(1,1)+x(1,n_cor_max+2))/2) )
    #         intens2 = I( floor((y(1,lengthl)+y(1,lengthl-n_cor_max-1))/2), floor((x(1,lengthl)+x(1,lengthl-n_cor_max-2))/2) )
    #         I(10,10)
    #         if intens1 > intens2
    #             #We need to reverse the numbering
    #             xTemp = x;
    #             yTemp = y;
    #             for i = 1:1:lengthl
    #                 x[i] = xTemp(lengthl+1 - i);
    #                 y[i] = yTemp(lengthl+1 - i);
    #             end
    #         end
    #     end


def _fix_numbering_direction(x: np.ndarray, y: np.ndarray, n_sq_x: int, n_sq_y: int) -> (np.ndarray, np.ndarray):
    ###########################################################################
    # Check, whether the numbering increases along the longer or the shorter pattern dimension:
    ###########################################################################
    n_cor_min = min(n_sq_x + 1, n_sq_y + 1)
    n_cor_max = max(n_sq_x + 1, n_sq_y + 1)
    n = n_cor_max * n_cor_min - 1
    hmin = np.zeros(n)
    hmin[:n:n_cor_min] = 1
    hmax = np.zeros(n)
    hmax[:n:n_cor_max] = 1
    dxy = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    pmin = np.convolve(hmin, dxy).max()
    pmax = np.convolve(hmax, dxy).max()
    inc_dir = n_cor_min if pmin > pmax else n_cor_max
    # Rearrange numbering from starting point
    xTemp = x.copy()
    yTemp = y.copy()
    n = len(x)
    area = np.empty(n - 1)
    for i in range(n - 1):
        xborder = np.hstack([
            xTemp[:inc_dir - 1],
            xTemp[inc_dir - 1: -inc_dir:inc_dir],
            xTemp[-inc_dir + 1:][::-1],
            xTemp[:-inc_dir + 1:inc_dir][::-1]
        ])
        yborder = np.hstack([
            yTemp[:inc_dir - 1],
            yTemp[inc_dir - 1:-inc_dir:inc_dir],
            yTemp[-inc_dir + 1:][::-1],
            yTemp[:-inc_dir + 1:inc_dir][::-1]
        ])
        area[i] = np.abs(np.trapz(xborder, yborder))
        xTemp = np.roll(xTemp, -1)
        yTemp = np.roll(yTemp, -1)
    shift = np.argmax(area)
    x = np.roll(x, shift)
    y = np.roll(y, shift)
    return x, y


def _fix_nonsquare_board_ambiguity(x: np.ndarray, y: np.ndarray, n_sq_x: int, n_sq_y: int) -> (np.ndarray, np.ndarray):
    # This algorithm was first designed to be used with square patterns only,
    # where the starting corner is meaningless, since every orientation +n*90
    # degrees is structurally equivalent.
    # With the extension to non-square checker boards, this is no longer the
    # case and, depending on the specific board, 2 or 4 orientations can be
    # distinguished.

    # Check, whether the numbering increases along the longer pattern dimension:
    # n_cor_min = min(n_sq_x+1, n_sq_y+1);
    # n_cor_max = max(n_sq_x+1, n_sq_y+1);
    # dist1 = (x(1,n_cor_min)-x(1,n_cor_min+1))^2 + (y(1,n_cor_min)-y(1,n_cor_min+1))^2;
    # dist2 = (x(1,n_cor_max)-x(1,n_cor_max+1))^2 + (y(1,n_cor_max)-y(1,n_cor_max+1))^2;
    n_cor_x = n_sq_x + 1
    n_cor_y = n_sq_y + 1
    dist1 = (x[n_cor_x - 1] - x[n_cor_x])**2 + (y[n_cor_x - 1] - y[n_cor_x])**2
    dist2 = (x[n_cor_y - 1] - x[n_cor_y])**2 + (y[n_cor_y - 1] - y[n_cor_y])**2
    if dist1 > dist2:
        # We have it wrongly numbered, renumber
        x_temp = x.copy()
        y_temp = y.copy()
        iter_offset = 0
        for i in range(len(x)):
            j = i % n_cor_y + 1
            x[i] = x_temp[j * n_cor_x - iter_offset - 1]
            y[i] = y_temp[j * n_cor_x - iter_offset - 1]
            if j * n_cor_x > n_cor_x * (n_cor_y - 1):
                iter_offset += 1
    return x, y


def _calculate_bounds(x: np.ndarray, y: np.ndarray, img_shape):
    # Calculate board bounds on image for visualization
    mask = (x >= 0) & (y >= 0)
    max_x = x[mask].max()
    min_x = x[mask].min()
    max_y = y[mask].max()
    min_y = y[mask].min()
    num_found_corners = mask.sum()
    num_corners = len(x)
    min_x -= (max_x - min_x) * (1 - num_found_corners / num_corners + 0.2)
    max_x += (max_x - min_x) * (1 - num_found_corners / num_corners + 0.2)
    min_y -= (max_y - min_y) * (1 - num_found_corners / num_corners + 0.2)
    max_y += (max_y - min_y) * (1 - num_found_corners / num_corners + 0.2)
    h, w = img_shape
    min_x = max(min_x, 0)
    max_x = min(max_x, w)
    min_y = max(min_y, 0)
    max_y = min(max_y, h)
    return min_x, max_x, min_y, max_y


def visualize_corners(x, y, img, title, n_sq_y, corner_color='red', axis_color='lime', cmap='gray'):
    img = np.array(img)
    fig, ax = plt.subplots()
    ax.set_title(title, dict(fontweight='bold'))
    ax.imshow(img, cmap=cmap, interpolation='bilinear')
    draw_axes(ax, y, x, n_sq_y, axis_color)
    ax.plot(x, y, '+', c=corner_color, linewidth=2)
    for i in range(len(x)):
        ax.text(x[i] + 3, y[i] + 3, str(i + 1), dict(color=corner_color))
    min_x, max_x, min_y, max_y = _calculate_bounds(x, y, img.shape)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(max_y, min_y)
    plt.tight_layout()
    plt.show()
