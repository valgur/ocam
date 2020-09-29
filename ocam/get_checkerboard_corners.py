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
# 04.03.2014 by Steffen Urban
# this is a modified file from
# Davide Scaramuzzas Toolbox OcamCalib
# filename: get_checkerboard_corners.m
# Code was changed to force subpixel corners

import subprocess
from glob import glob

from .calib_data import CalibData
from .cornerfinder import cornerfinder
from .libsmop import *


def call_find_corners(img_file, n_sq_x, n_sq_y):
    img_file = os.path.abspath(img_file)
    cur_dir = os.curdir
    script_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(script_dir, 'autoCornerFinder'))
    try:
        for f in glob('cToMatlab/*.txt'):
            os.remove(f)

        # Tell the automatic corner extractor, which image file to process
        with open('pictures.txt', 'w') as f:
            f.write(img_file)

        # Call the automatic corner extraction algorithm
        # Width and height in the algorithm are defined differently than
        # in this toolbox: Its the number of inernal corners instead of 
        # internal quadrangles.
        # -m specifies the number of corners the algorithm has to find, before it
        # terminates
        # If callback = -1   -> An error occured, automatic corner finding is aborted
        #            =  0   -> Not enough corners have been found, add some manually
        #            =  1   -> Enough corners have been found for calibration

        retcode = subprocess.call([
            './FindCorners',
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
        cornerInfo = np.loadtxt('cToMatlab/cornerInfo.txt')
        cornersX = np.loadtxt('cToMatlab/cornersX.txt')
        cornersY = np.loadtxt('cToMatlab/cornersY.txt')
    finally:
        for f in glob('cToMatlab/*.txt'):
            os.remove(f)
        os.chdir(cur_dir)
    return cornerInfo, cornersX, cornersY


@function
def get_checkerboard_corners(kk: int, use_corner_find: bool, calib_data: CalibData):
    # INITIALIZATIONS
    ###########################################################################
    print(f'Processing image {calib_data.imgs[kk]}...')
    img = calib_data.read_image(kk)

    cornerInfo, cornersX, cornersY = call_find_corners(
        os.path.abspath(calib_data.imgs[kk]), calib_data.n_sq_x, calib_data.n_sq_y)
    cornerInfo = matlabarray(cornerInfo)
    cornersX = matlabarray(cornersX)
    cornersY = matlabarray(cornersY)

    # VARIABLE DEFINITIONS
    ###########################################################################
    numCorners = (calib_data.n_sq_x + 1) * (calib_data.n_sq_y + 1)
    numOfFoundCorners = 0
    cornerNumber = 0
    startingCorner = copy([])
    deltaCols = 0
    startingCornerSet = false
    # Index of the corner with the smallest index
    min_i = copy([])
    min_j = copy([])
    # Used for image zoom in
    size1, size2 = img.shape
    min_x = max(size1, size2) + 1
    max_x = 0
    min_y = max(size1, size2) + 1
    max_y = 0

    # INPUT CORNER VALUE AND MATRIX SIZE ADAPTATIONS
    ###########################################################################
    # If there are entire rows of zeros at the end of the file, this means that
    # the cornerfinder is not entirely sure about the location of the board.
    # Then eventually additional corners are needed in the front.
    iSave, _ = find(cornersX >= 0, nargout=2)
    iSave = max(iSave)
    # Add at least one row! This is because the columns are maybe ambiguous, 
    # and then we don't know whether they belong to the current row or to the
    # next one...
    if cornerInfo[0] - iSave >= 0:
        cornersX = concat([
            [-1 * ones(cornerInfo[0] - iSave + 1, cornerInfo[1])],
            [cornersX]
        ])
        cornersY = concat([
            [-1 * ones(cornerInfo[0] - iSave + 1, cornerInfo[1])],
            [cornersY]
        ])

    rows, cols = cornersX.shape
    # Add one pixel to every non "-1" value, since Matlab starts numbering
    # at one, whereas c++ starts at zero.
    flagStart = True
    for i in arange(1, rows).flat:
        for j in arange(1, cols).flat:
            # Define the starting corner as the first found corner
            # which has a neighbor corner which was also found
            if j != cols and flagStart == true:
                if cornersX(i, j) >= 0 and cornersX(i, j + 1) >= 0:
                    startingCorner = concat([i, j])
                    flagStart = false
            if cornersX(i, j) >= 0:
                cornersX[i, j] = cornersX(i, j) + 1
                numOfFoundCorners += 1
                # found corners. ->Needed further down
                if cornersX(i, j) > max_x:
                    max_x = cornersX(i, j)
                if cornersX(i, j) < min_x:
                    min_x = cornersX(i, j)
            if cornersY(i, j) >= 0:
                cornersY[i, j] = cornersY(i, j) + 1
                if cornersY(i, j) > max_y:
                    max_y = cornersY(i, j)
                if cornersY(i, j) < min_y:
                    min_y = cornersY(i, j)

    # PREPARATIONS FOR PROPER PLOT ZOOM-IN
    ###########################################################################
    min_x -= (max_x - min_x) * (1 - numOfFoundCorners / numCorners + 0.2)
    max_x += (max_x - min_x) * (1 - numOfFoundCorners / numCorners + 0.2)
    min_y -= (max_y - min_y) * (1 - numOfFoundCorners / numCorners + 0.2)
    max_y += (max_y - min_y) * (1 - numOfFoundCorners / numCorners + 0.2)
    min_x = max(min_x, 0)
    max_x = min(max_x, size2)
    min_y = max(min_y, 0)
    max_y = min(max_y, size1)

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
    #     set(h, 'FontWeight', 'bold')
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
    #     set(findobj('type', 'text'), 'color', 'red'); 
    #     hold off;
    #       
    #     #Get user input on behalf of the board orientation
    #     numCornersDirZeroOne = input(['The automatic corner finder was not able to decide how the pattern is oriented. Please indicate the number of corners in direction of corners [0 -> 1]: ']);
    #     #We look in row direction
    #     deltaCols = cornerInfo[1] - numCornersDirZeroOne;                 
    # end

    # DRAW AND NUMBER FOUND CORNERS
    ###########################################################################
    # Draw the found corners onto the image and number them 
    # Define the first encountered found corner as number "zero"
    PlotXxiX = copy([])
    PlotXxiY = copy([])
    PlotCornersX = copy([])
    PlotCornersY = copy([])
    PlotCornerNumber = copy([])
    for i in arange(1, rows).flat:
        for j in arange(1, (cols - deltaCols)).flat:
            if cornersX(i, j) != - 1:
                # Save for plotting later
                PlotCornersX = concat([PlotCornersX, cornersX(i, j)])
                PlotCornersY = concat([PlotCornersY, cornersY(i, j)])
                PlotCornerNumber = concat([PlotCornerNumber, cornerNumber])
                if use_corner_find:
                    # Apply the corner finder
                    #   [xxi] = cornerfinder([cornersX(i,j);cornersY(i,j)],I,winty,wintx);
                    ## ================ 
                    #  added(changed code steffen urban
                    wintx = ceil(calib_data.ocam_model.height / 100)
                    winty = ceil(calib_data.ocam_model.height / 100)
                    xxi = cornerfinder(concat([[cornersX(i, j)], [cornersY(i, j)]]), img, winty, wintx)
                    cornersX[i, j] = xxi[1]
                    cornersY[i, j] = xxi[2]
                    # Save for plotting later
                    PlotXxiX = concat([PlotXxiX, xxi[1]])
                    PlotXxiY = concat([PlotXxiY, xxi[2]])
                if cornerNumber == 0:
                    # Change starting corner to this one
                    # Needed further down
                    startingCorner = concat([i, j])
                    startingCornerSet = true
            if startingCornerSet == true:
                cornerNumber += 1

    # If no negative corners were designated, then "min_i" and "min_j"
    # are still empty
    if isempty(min_i):
        min_i = startingCorner[1]
        min_j = startingCorner[2]

    # SAVE CORNER INFORMATION IN 2 VECTORS
    ###########################################################################
    # Save all corners in two arrays for further processing
    # by other functions
    x = copy([])
    y = copy([])
    # Start with the smallest found corner and then append the larger ones
    iteration = 0
    while True:
        # Iterators
        i = min_i + floor(iteration / (cols - deltaCols))
        j = mod((min_j - 1 + iteration), cols - deltaCols) + 1
        iteration += 1
        x = concat([x, cornersX(i, j)])
        y = concat([y, cornersY(i, j)])
        if iteration >= numCorners or i > rows:
            break

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
    #                 x(i) = xTemp(lengthl+1 - i);
    #                 y(i) = yTemp(lengthl+1 - i);
    #             end
    #         end
    #     end

    # REPOSITION CORNERS (IF NEEDED)
    ###########################################################################
    # Allow the user (if needed) to reposition any of the corners
    sizex = x.shape[1]
    iNumber = arange(1, sizex)
    ###########################################################################
    # Check, whether the numbering increases along the longer or the shorter pattern dimention:
    ###########################################################################
    n_cor_min = min(calib_data.n_sq_x + 1, calib_data.n_sq_y + 1)
    n_cor_max = max(calib_data.n_sq_x + 1, calib_data.n_sq_y + 1)
    hmin = zeros(1, n_cor_max * n_cor_min - 1)
    hmin[arange(1, end(), n_cor_min)] = 1
    hmax = zeros(1, n_cor_max * n_cor_min - 1)
    hmax[arange(1, end(), n_cor_max)] = 1
    dxy = sqrt(diff(x)**2 + diff(y)**2)
    pmin = conv(hmin, dxy)
    pmin = max(ravel(pmin))
    pmax = conv(hmax, dxy)
    pmax = max(ravel(pmax))
    if pmin[1] > pmax[1]:
        inc_dir = n_cor_min
        oth_dir = n_cor_max
    else:
        inc_dir = n_cor_max
        oth_dir = n_cor_min

    # Rearrange numbering from starting point
    xTemp = x
    yTemp = y
    area = copy([])
    for i in arange(1, length(x) - 1).flat:
        xborder = concat([
            xTemp(1, arange(1, inc_dir - 1)),
            xTemp(1, arange(inc_dir, end() - inc_dir, inc_dir)),
            xTemp(1, arange(end(), end() - inc_dir + 2, - 1)),
            xTemp(1, arange(end() - inc_dir + 1, 1, - inc_dir))
        ])
        yborder = concat([
            yTemp(1, arange(1, inc_dir - 1)),
            yTemp(1, arange(inc_dir, end() - inc_dir, inc_dir)),
            yTemp(1, arange(end(), end() - inc_dir + 2, - 1)),
            yTemp(1, arange(end() - inc_dir + 1, 1, - inc_dir))
        ])
        area[i] = abs(trapz(xborder, yborder))
        xTemp = concat([xTemp(1, arange(2, end())), xTemp(1, 1)])
        yTemp = concat([yTemp(1, arange(2, end())), yTemp(1, 1)])

    shift = find(area == max(area)) - 1
    if shift > 0:
        x = concat([x(1, arange(1 + shift, end())), x(1, arange(1, shift))])
        y = concat([y(1, arange(1 + shift, end())), y(1, arange(1, shift))])

    ###########################################################################
    # ASSIGN STARTING CORNER
    ###########################################################################
    # This algorithm was first designed to be used with square patterns only,
    # where the starting corner is meaningless, since every orientation +n*90
    # degrees is structurally equivalent.
    # With the extension to non-square checker boards, this is no longer the
    # case and, depending on the specific board, 2 or 4 orientations can be
    # distinguished.
    # 1. Check, whether we are dealing with a non-square board:
    if calib_data.n_sq_x != calib_data.n_sq_y:
        # 2. Check, whether the numbering increases along the longer pattern
        # dimention:
        #     n_cor_min = min(n_sq_x+1, n_sq_y+1);
        #     n_cor_max = max(n_sq_x+1, n_sq_y+1);
        #     dist1 = (x(1,n_cor_min)-x(1,n_cor_min+1))^2 + (y(1,n_cor_min)-y(1,n_cor_min+1))^2;
        #     dist2 = (x(1,n_cor_max)-x(1,n_cor_max+1))^2 + (y(1,n_cor_max)-y(1,n_cor_max+1))^2;
        n_cor_x = calib_data.n_sq_x + 1
        n_cor_y = calib_data.n_sq_y + 1
        dist1 = (x(1, n_cor_x) - x(1, n_cor_x + 1))**2 + (y(1, n_cor_x) - y(1, n_cor_x + 1))**2
        dist2 = (x(1, n_cor_y) - x(1, n_cor_y + 1))**2 + (y(1, n_cor_y) - y(1, n_cor_y + 1))**2
        if dist1 > dist2:
            # We have it wrongly numbered, renumber
            xTemp = x
            yTemp = y
            lengthl = x.shape[1]
            iterMult = n_cor_x
            iterOffset = 0
            for i in arange(1, lengthl).flat:
                j = mod(i - 1, n_cor_y) + 1
                x[i] = xTemp(j * iterMult - iterOffset)
                y[i] = yTemp(j * iterMult - iterOffset)
                if j * iterMult > n_cor_x * (n_cor_y - 1):
                    iterOffset += 1

    Yp_abs = x.T
    Xp_abs = y.T
    ###########################################################################
    # Visualize
    ###########################################################################
    # if 1:
    #     figure[2]
    #     imagesc(I)
    #     colormap(gray)
    #     title(concat(['Image ', str(kk)]))
    #     h = get(gca, 'title')
    #     set(h, 'FontWeight', 'bold')
    #     hold('on')
    #     plot(x, y, '+', 'color', 'red', 'linewidth', 2)
    #     text(x.T + 3, y.T + 3, str(iNumber.T))
    #     set(findobj('type', 'text'), 'color', 'red')
    #     axis(concat([min_x, max_x, min_y, max_y]))
    #     draw_axes(Xp_abs, Yp_abs, calib_data.n_sq_y)
    #     drawnow
    #     hold('off')
    #     #     close[2];

    ###########################################################################
    # Display finished message
    ###########################################################################
    print('Done')

    return Xp_abs, Yp_abs
