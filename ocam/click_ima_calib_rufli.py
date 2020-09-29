###########################################################################  
#   Copyright (C) 2007 MARTIN RUFLI
#   
#   Initially written by Martin Rufli and modified by Davide Scaramuzza
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
###########################################################################

from .calib_data import CalibData
from .libsmop import *


@function
def click_ima_calib_rufli(kk, use_corner_find, calib_data: CalibData):
    if exist('autoCornerFinder/cToMatlab/cornerInfo.txt', 'file'):
        os.remove('autoCornerFinder/cToMatlab/cornerInfo.txt')

    if exist('autoCornerFinder/cToMatlab/cornersX.txt', 'file'):
        os.remove('autoCornerFinder/cToMatlab/cornersX.txt')

    if exist('autoCornerFinder/cToMatlab/cornersY.txt', 'file'):
        os.remove('autoCornerFinder/cToMatlab/cornersY.txt')

    if exist('autoCornerFinder/cToMatlab/error.txt', 'file'):
        os.remove('autoCornerFinder/cToMatlab/error.txt')

    # INITIALIZATIONS
    ###########################################################################
    print('\nProcessing image %d...\n', kk, end='')
    I = calib_data.I[kk]
    # if ~(size(calib_data.wintx_,2)<kk),
    #     
    #     wintxkk = calib_data.wintx_{kk};
    #     
    #     if ~isempty(wintxkk) & ~isnan(wintxkk),
    #         
    #         calib_data.wintx = calib_data.wintx_{kk};
    #         calib_data.winty = calib_data.winty_{kk};
    #         
    #     end;
    # end;

    if use_corner_find:
        print('Using (wintx,winty)=(%d,%d) - Window size = %dx%d\n', calib_data.wintx, calib_data.winty,
              dot(2, calib_data.wintx) + 1, dot(2, calib_data.winty) + 1, end='')

    # EXTERNAL CORNER FINDING ALGORITHM CALL
    ###########################################################################
    # kk jumps over empty image slots (i.e. kk = 1->2->4->5, if image 3 was not 
    # loaded), whereas l(j,1).name does not. Therefore account for it!
    # if kk == indices( ima_numbers[1] );
    #     iter_r = ima_numbers[1];
    #     iter_rr= 2;
    # else
    #     iter_r = ima_numbers(iter_rr);
    #     iter_rr += 1;
    # end

    # Tell the automatic corner extractor, which image file to process
    fid = open('./autoCornerFinder/pictures.txt', 'w')
    fid.write('../%s' % (calib_data.imgs[kk]))
    fid.close()
    # iter_r += 1;    #!!!iter_r is defined in "click_calib.m"!!!

    # Call the automatic corner extraction algorithm
    # Width and height in the algorithm are defined differently than
    # in this toolbox: Its the number of inernal corners instead of 
    # internal quadrangles.
    # -m specifies the number of corners the algorithm has to find, before it
    # terminates
    # If callback = -1   ->An error occured, automatic corner finding is aborted
    #            =  0   ->Not enough corners have been found, add some manually
    #            =  1   ->Enough corners have been found for calibration

    # Visualization turned OFF
    os.chdir('autoCornerFinder')
    callString = (concat(['FindCorners -w ', str(calib_data.n_sq_x + 1), ' -h ', str(calib_data.n_sq_y + 1), ' -m ',
                          str(dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1))), ' pictures.txt']))
    if logical_not(ispc):
        callString = concat(['./', callString])

    # Visualization turned ON
    # callString = (['cd autoCornerFinder & FindCornersVisual -w ' str(n_sq_x+1) ' -h ' str(n_sq_y+1) ' -m ' str((n_sq_x+1) * (n_sq_y+1)) ' pictures.txt']);

    # Visualization turned ON and Saving of the images turned ON
    # WARNING: Does somehow not work under Windows 2000...
    # callString = (['cd autoCornerFinder & FindCornersVisualSave -w ' str(n_sq_x+1) ' -h ' str(n_sq_y+1) ' -m ' str((n_sq_x+1) * (n_sq_y+1)) ' pictures.txt']);

    # system('who');
    callBack = system(callString)
    os.chdir('..')
    # Do error checking
    if callBack == - 1:
        # Display the error message
        disp('During corner finding an error occured:')
        filename = 'autoCornerFinder/cToMatlab/error.txt'
        fid = open(filename, 'r')
        line = fgetl(fid)
        fid.close()
        disp(line)
        disp('Please restart "Extract grid corners" or remove this image from the input dataset.')
        return

    # Open the corner size information file
    filename = 'autoCornerFinder/cToMatlab/cornerInfo.txt'
    fid = open(filename, 'r')
    cornerInfo = fscanf(fid, '%g %g', concat([1, 2]))
    fid.close()
    # Open the files with the found corners
    filename = 'autoCornerFinder/cToMatlab/cornersX.txt'
    fid = open(filename, 'r')
    cornersX = fscanf(fid, '%g %g', concat([cornerInfo[2], cornerInfo[1]]))
    cornersX = cornersX.T
    fid.close()
    filename = 'autoCornerFinder/cToMatlab/cornersY.txt'
    fid = open(filename, 'r')
    cornersY = fscanf(fid, '%g %g', concat([cornerInfo[2], cornerInfo[1]]))
    cornersY = cornersY.T
    fid.close()
    # VARIABLE DEFINITIONS
    ###########################################################################
    numCorners = dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1))
    numCornerThreshold = dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1))
    numOfFoundCorners = 0
    cornerNumber = 0
    cornersNumbering = dot(- 1, ones(cornerInfo[1], cornerInfo[2]))
    startingCorner = copy([])
    nextFoundCorner = concat([- 1, - 1])
    deltaCols = 0
    startingCornerSet = false
    # Index of the corner with the smallest index
    min_i = copy([])
    min_j = copy([])
    # Used for image zoom in
    size1, size2 = size(I, nargout=2)
    min_x = max(size1, size2) + 1
    max_x = 0
    min_y = max(size1, size2) + 1
    max_y = 0
    # INPUT CORNER VALUE AND MATRIX SIZE ADAPTATIONS
    ###########################################################################
    # If there are entire rows of zeros at the end of the file, this means that
    # the cornerfinder is not entirely shure about the location of the board.
    # Then eventually additional corners are needed in the front.
    iSave, dummy = find(cornersX >= 0, nargout=2)
    iSave = max(iSave)
    # Add at least one row! This is because the columns are maybe ambiguous, 
    # and then we don't know whether they belong to the current row or to the
    # next one...
    if (cornerInfo[1] - iSave >= 0):
        cornersX = concat([[dot(- 1, ones(cornerInfo[1] - iSave + 1, cornerInfo[2]))], [cornersX]])
        cornersY = concat([[dot(- 1, ones(cornerInfo[1] - iSave + 1, cornerInfo[2]))], [cornersY]])

    rows, cols = size(cornersX, nargout=2)
    # Add one pixel to every non "-1" value, since Matlab starts numbering
    # at one, whereas c++ starts at zero.
    flagStart = true
    for i in arange(1, rows, 1).flat:
        for j in arange(1, cols, 1).flat:
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
    min_x -= dot((max_x - min_x), (1 - numOfFoundCorners / numCorners + 0.2))
    max_x += dot((max_x - min_x), (1 - numOfFoundCorners / numCorners + 0.2))
    min_y -= dot((max_y - min_y), (1 - numOfFoundCorners / numCorners + 0.2))
    max_y += dot((max_y - min_y), (1 - numOfFoundCorners / numCorners + 0.2))
    min_x = max(min_x, 0)
    max_x = min(max_x, size2)
    min_y = max(min_y, 0)
    max_y = min(max_y, size1)
    # ORIENTATION AMBIGUITY?
    ###########################################################################
    # Decide whether the position and orientation of the checkerboard
    # can be unambiguously determined
    if min((calib_data.n_sq_y + 1), (calib_data.n_sq_x + 1)) != min(cornerInfo[1], cornerInfo[2]):
        if flagStart == true:
            # DISPLAY AN ERROR AND RETURN
            pass
        # There is some ambiguity, get rid of it
        figure[2]
        image(I)
        colormap(calib_data.map)
        set(2, 'color', concat([1, 1, 1]))
        title(cellarray([concat(['Image ', str(kk)])]))
        h = get(gca, 'title')
        set(h, 'FontWeight', 'bold')
        axis(concat([min_x, max_x, min_y, max_y]))
        figure[2]
        hold('on')
        i = startingCorner[1]
        j = startingCorner[2]
        # Plot the starting corner and its neighbor
        plot(cornersX(i, j), cornersY(i, j), '+', 'color', 'red', 'linewidth', 2)
        plot(cornersX(i, j + 1), cornersY(i, j + 1), '+', 'color', 'red', 'linewidth', 2)
        text(cornersX(i, j) + 3, cornersY(i, j) + 3, str[0])
        text(cornersX(i, j + 1) + 3, cornersY(i, j + 1) + 3, str[1])
        set(findobj('type', 'text'), 'color', 'red')
        hold('off')
        numCornersDirZeroOne = input(concat([
            'The automatic corner finder was not able to decide how the pattern is oriented. Please indicate the number of corners in direction of corners [0 -> 1]: ']))
        deltaCols = cornerInfo[2] - numCornersDirZeroOne

    # DRAW AND NUMBER FOUND CORNERS
    ###########################################################################
    # Draw the found corners onto the image and number them 
    # Define the first encountered found corner as number "zero"
    PlotXxiX = copy([])
    PlotXxiY = copy([])
    PlotCornersX = copy([])
    PlotCornersY = copy([])
    PlotCornerNumber = copy([])
    for i in arange(1, rows, 1).flat:
        for j in arange(1, (cols - deltaCols), 1).flat:
            if cornersX(i, j) != - 1:
                # Save for plotting later
                PlotCornersX = concat([PlotCornersX, cornersX(i, j)])
                PlotCornersY = concat([PlotCornersY, cornersY(i, j)])
                PlotCornerNumber = concat([PlotCornerNumber, cornerNumber])
                if use_corner_find:
                    # Apply the corner finder
                    xxi = cornerfinder(concat([[cornersX(i, j)], [cornersY(i, j)]]), I, winty, wintx)
                    PlotXxiX = concat([PlotXxiX, xxi[1]])
                    PlotXxiY = concat([PlotXxiY, xxi[2]])
                if cornerNumber == 0:
                    # Change starting corner to this one
                    # Needed further down
                    startingCorner = concat([i, j])
                    startingCornerSet = true
            if startingCornerSet == true:
                cornerNumber += 1

    figure[2]
    image(I)
    colormap(calib_data.map)
    set(2, 'color', concat([1, 1, 1]))
    title(cellarray([[concat(['Image ', str(kk)])],
                     [concat([str(numOfFoundCorners), ' / ', str(numCorners), ' corner have been found.'])],
                     [concat(['Press ENTER to continue.'])]]))
    h = get(gca, 'title')
    set(h, 'FontWeight', 'bold')
    axis(concat([min_x, max_x, min_y, max_y]))
    # Plot the original corners
    figure[2]
    hold('on')
    plot(PlotCornersX, PlotCornersY, '+', 'color', 'red', 'linewidth', 2)
    if use_corner_find:
        # Plot the "corner finder enhanced" corners
        plot(PlotXxiX, PlotXxiY, '+', 'color', 'blue', 'linewidth', 2)

    text(PlotCornersX.T + 3, PlotCornersY.T + 3, str(PlotCornerNumber.T))
    set(findobj('type', 'text'), 'color', 'red')
    hold('off')
    pause
    # ADD NEW CORNERS
    ###########################################################################
    # Only do this, if we still need to add some corners
    if numCorners != numOfFoundCorners:
        disp('Press ENTER and then click on the corner whose number is highlighted in the title bar.')
        disp(
            'Corner selection starts in increasing order. For changing mode between increasing and decreasing order, right click on the image.')
        pause
        figure[2]
        hold('on')
        mode = 1
        nextCorner = 0
        cornerNumberMin = 0
        iteration = startingCorner[2]
        while true:

            # Iterators & update (i,j)
            if mode == 1:
                i = startingCorner[1] + floor((iteration - 1) / (cols - deltaCols))
                j = mod((iteration - 1), cols - deltaCols) + 1
                iteration += 1
            else:
                i = startingCorner[1] + floor((iteration - 1) / (cols - deltaCols))
                j = mod((iteration - 1), cols - deltaCols) + 1
                iteration -= 1
            # Check whether i or j are out of bounds
            # If yes, switch mode
            if (i <= logical_or(0, j) <= logical_or(0, i) > logical_or(rows, j) > cols - deltaCols):
                mode = dot(mode, - 1)
                i = startingCorner[1]
                j = startingCorner[2]
                # iteration = 0;
                iteration = startingCorner[2]
                nextCorner = 0
                continue
            # Continue, if corner is already labeled
            if cornersX(i, j) != - 1:
                if mode == 1:
                    nextCorner += 1
                else:
                    nextCorner -= 1
                continue
            figure[2]
            if (numCornerThreshold - numOfFoundCorners) > 1:
                title(cellarray([[concat(['Image ', str(kk)])],
                                 [concat([str(numCornerThreshold - numOfFoundCorners), ' corner are missing.'])],
                                 [concat(['Please place corner no. ', str(nextCorner), ' on the plot.'])]]))
            else:
                title(cellarray([[concat(['Image ', str(kk)])],
                                 [concat([str(numCornerThreshold - numOfFoundCorners), ' corner is missing.'])],
                                 [concat(['Please place corner no. ', str(nextCorner), ' on the plot.'])]]))
            h = get(gca, 'title')
            set(h, 'FontWeight', 'bold')
            xi, yi, button = ginput3(1, nargout=3)
            if button > 1:
                mode = dot(mode, - 1)
                i = startingCorner[1]
                j = startingCorner[2]
                # iteration = 0;
                nextCorner = 0
                iteration = startingCorner[2]
                continue
            if use_corner_find:
                # Use corner enhancer
                xxi = cornerfinder(concat([[xi], [yi]]), I, winty, wintx)
                xi = xxi[1]
                yi = xxi[2]
            figure[2]
            plot(xi, yi, '+', 'color', 'red', 'linewidth', 2)
            text(xi + 3, yi + 3, str(nextCorner))
            set(findobj('type', 'text'), 'color', 'red')
            cornersX[i, j] = xi
            cornersY[i, j] = yi
            if nextCorner <= cornerNumberMin:
                min_i = i
                min_j = j
            # Adjust "nextCorner", depending on the current mode
            if mode == 1:
                nextCorner += 1
            else:
                nextCorner -= 1
            # If we are here, a new corner has been found
            # Increase the found corner count
            numOfFoundCorners += 1
            if numOfFoundCorners >= numCornerThreshold:
                break

        hold('off')

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
    while true:

        # Iterators
        i = min_i + floor(iteration / (cols - deltaCols))
        j = mod((min_j - 1 + iteration), cols - deltaCols) + 1
        iteration += 1
        x = concat([x, cornersX(i, j)])
        y = concat([y, cornersY(i, j)])
        if (iteration >= logical_or(numCorners, i) > rows):
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
    dummy, sizex = size(x, nargout=2)
    iNumber = arange(1, sizex, 1)
    reposition_corner = input('Would you like to reposition any of the assigned corners ([] = yes, other = no)?', 's')
    if isempty(reposition_corner):
        figure[2]
        image(I)
        colormap(calib_data.map)
        set(2, 'color', concat([1, 1, 1]))
        title(cellarray([[concat(['Image ', str(kk)])], [concat(['Next you can reposition badly placed corners.'])],
                         [concat(['Press ENTER to continue.'])]]))
        h = get(gca, 'title')
        set(h, 'FontWeight', 'bold')
        axis(concat([min_x, max_x, min_y, max_y]))
        figure[2]
        hold('on')
        dummy, sizex = size(x, nargout=2)
        iNumber = arange(1, sizex, 1)
        while true:

            # Display the (updated) corners
            figure[2]
            plot(x, y, '+', 'color', 'red', 'linewidth', 2)
            text(x.T + 3, y.T + 3, str(iNumber.T))
            set(findobj('type', 'text'), 'color', 'red')
            figure[2]
            title(cellarray([[concat(['Image ', str(kk)])], [concat(['Left click on a corner to replace it.'])],
                             [concat(['Right click anywhere to quit replacement mode.'])]]))
            h = get(gca, 'title')
            set(h, 'FontWeight', 'bold')
            xdel, ydel, button = ginput3(1, nargout=3)
            if button > 1:
                break
            xdel = dot(xdel, ones(dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1)), 1))
            ydel = dot(ydel, ones(dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1)), 1))
            distMatrix = concat([xdel, ydel]) - concat([x.T, y.T])
            distMatrix = dot(distMatrix.T, distMatrix.T)
            nearestCornerID = find(floor(distMatrix.T) == min(floor(distMatrix.T)))
            if (nearestCornerID == dot((calib_data.n_sq_x + 1), (calib_data.n_sq_y + 1))):
                x[nearestCornerID] = x(nearestCornerID - 1)
                y[nearestCornerID] = y(nearestCornerID - 1)
                iNumber[nearestCornerID] = iNumber(nearestCornerID - 1)
            else:
                x[nearestCornerID] = x(nearestCornerID + 1)
                y[nearestCornerID] = y(nearestCornerID + 1)
                iNumber[nearestCornerID] = iNumber(nearestCornerID + 1)
            figure[2]
            image(I)
            colormap(calib_data.map)
            set(2, 'color', concat([1, 1, 1]))
            title(cellarray([[concat(['Image ', str(kk)])], [concat(['Left click on the desired new location.'])],
                             [concat([' '])]]))
            h = get(gca, 'title')
            set(h, 'FontWeight', 'bold')
            figure[2]
            hold('on')
            plot(x, y, '+', 'color', 'red', 'linewidth', 2)
            text(x.T + 3, y.T + 3, str(iNumber.T))
            set(findobj('type', 'text'), 'color', 'red')
            xnew, ynew, button = ginput3(1, nargout=3)
            x[nearestCornerID] = xnew
            y[nearestCornerID] = ynew
            iNumber = arange(1, sizex, 1)

        hold('off')

    ###########################################################################
    # Check, whether the numbering increases along the longer or the shorter pattern dimention:
    ###########################################################################
    n_cor_min = min(calib_data.n_sq_x + 1, calib_data.n_sq_y + 1)
    n_cor_max = max(calib_data.n_sq_x + 1, calib_data.n_sq_y + 1)
    hmin = zeros(1, dot(n_cor_max, n_cor_min) - 1)
    hmin[arange(1, end(), n_cor_min)] = 1
    hmax = zeros(1, dot(n_cor_max, n_cor_min) - 1)
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
        xborder = concat([xTemp(1, arange(1, inc_dir - 1)), xTemp(1, arange(inc_dir, end() - inc_dir, inc_dir)),
                          xTemp(1, arange(end(), end() - inc_dir + 2, - 1)),
                          xTemp(1, arange(end() - inc_dir + 1, 1, - inc_dir))])
        yborder = concat([yTemp(1, arange(1, inc_dir - 1)), yTemp(1, arange(inc_dir, end() - inc_dir, inc_dir)),
                          yTemp(1, arange(end(), end() - inc_dir + 2, - 1)),
                          yTemp(1, arange(end() - inc_dir + 1, 1, - inc_dir))])
        area[i] = abs(trapz(xborder, yborder))
        xTemp = concat([xTemp(1, arange(2, end())), xTemp(1, 1)])
        yTemp = concat([yTemp(1, arange(2, end())), yTemp(1, 1)])

    shift = find(area == max(area)) - 1
    if shift > 0:
        x = concat([x(1, arange(1 + shift, end())), x(1, arange(1, shift))])
        y = concat([y(1, arange(1 + shift, end())), y(1, arange(1, shift))])

    ###########################################################################
    # Visualize
    ###########################################################################
    if 0:
        figure[2]
        image(I)
        colormap(calib_data.map)
        hold('on')
        plot(x, y, '+', 'color', 'red', 'linewidth', 2)
        text(x.T + 3, y.T + 3, str(iNumber.T))
        set(findobj('type', 'text'), 'color', 'red')
        axis(concat([min_x, max_x, min_y, max_y]))

    ###########################################################################

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
            dummy, lengthl = size(x, nargout=2)
            iterMult = n_cor_x
            iterOffset = 0
            for i in arange(1, lengthl, 1).flat:
                j = mod(i - 1, n_cor_y) + 1
                x[i] = xTemp(dot(j, iterMult) - iterOffset)
                y[i] = yTemp(dot(j, iterMult) - iterOffset)
                if dot(j, iterMult) > dot(n_cor_x, (n_cor_y - 1)):
                    iterOffset += 1

    calib_data.Yp_abs[:, :, kk] = x.T
    calib_data.Xp_abs[:, :, kk] = y.T
    ###########################################################################
    # Visualize
    ###########################################################################
    if 1:
        figure[2]
        image(I)
        colormap(calib_data.map)
        title(cellarray(
            [[concat(['Image ', str(kk)])], [concat(['The corners have been renumbered in the right order.'])],
             [concat(['Press ENTER to continue.'])]]))
        h = get(gca, 'title')
        set(h, 'FontWeight', 'bold')
        hold('on')
        plot(x, y, '+', 'color', 'red', 'linewidth', 2)
        text(x.T + 3, y.T + 3, str(iNumber.T))
        set(findobj('type', 'text'), 'color', 'red')
        axis(concat([min_x, max_x, min_y, max_y]))
        draw_axes(calib_data.Xp_abs[:, :, kk], calib_data.Yp_abs[:, :, kk],
                  calib_data.n_sq_y)
        print('Press ENTER to continue.')
        pause
        close_[2]

    ###########################################################################
    # os.remove FILES
    ###########################################################################
    # We os.remove the interface files between Matlab and c++, in order to prevent
    # reloading old data in case of some errors.
    os.remove('autoCornerFinder/cToMatlab/cornerInfo.txt')
    os.remove('autoCornerFinder/cToMatlab/cornersX.txt')
    os.remove('autoCornerFinder/cToMatlab/cornersY.txt')
    os.remove('autoCornerFinder/cToMatlab/error.txt')
    ###########################################################################
    # END OF CODE
    return


if __name__ == '__main__':
    pass
