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

from .calib_data import CalibData
from .libsmop import *
from .world2cam import world2cam


#################### REPROJECT ON THE IMAGES ########################
@function
def reproject_calib(calib_data: CalibData):
    varargin = reproject_calib.varargin
    nargin = reproject_calib.nargin

    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    # if isempty(calib_data.no_image),
    #   no_image = 0;
    # end;

    if logical_or(isempty(calib_data.width), isempty(calib_data.height)):
        print('WARNING: No image size (width,height) available. Setting width=640 and height=480')
        calib_data.width = 640
        calib_data.height = 480

    # Color code for each image:

    colors = 'brgkcm'
    # Reproject the patterns on the images, and compute the pixel errors:

    # Reload the images if necessary
    if calib_data.n_imgs != 0:
        if isempty(calib_data.ocam_model.ss):
            print('Need to calibrate before showing image reprojection. Maybe need to load Calib_Results.mat file.')
            return

    if calib_data.n_imgs != 0:
        if false:
            if isempty(calib_data.ind_active[1]) or size(calib_data.I) < calib_data.ind_active[1] or isempty(
                    calib_data.I[calib_data.ind_active[1]]):
                n_ima_save = calib_data.n_imgs
                active_images_save = calib_data.active_images
                ima_read_calib(calib_data)
                calib_data.n_imgs = n_ima_save
                calib_data.active_images = active_images_save
                if calib_data.no_image_file:
                    print('WARNING: Do not show the original images')
        else:
            calib_data.no_image_file = 1

    if (isempty(calib_data.ocam_model.c) or isempty(calib_data.ocam_model.d) or isempty(calib_data.ocam_model.e)):
        calib_data.ocam_model.c = 1
        calib_data.ocam_model.d = 0
        calib_data.ocam_model.e = 0

    for kk in calib_data.ima_proc.flat:
        if (size(calib_data.I, 1) >= kk) and logical_not(isempty(calib_data.I[kk])):
            I = calib_data.I[kk]
        else:
            I = dot(255, ones(calib_data.height, calib_data.width))
        xx = dot(calib_data.RRfin[:, :, kk],
                 concat([[calib_data.Xt.T], [calib_data.Yt.T], [ones(size(calib_data.Xt.T))]]))
        m = world2cam(xx, calib_data.ocam_model)
        xp = m[1, :]
        yp = m[2, :]
        figure(5 + kk)
        image(I)
        hold('on')
        colormap(gray[256])
        hold('on')
        title(concat(['Image ', str(kk), ' - Image points (+) and reprojected grid points (o)']))
        plot(calib_data.Yp_abs[:, :, kk], calib_data.Xp_abs[:, :, kk], 'r+')
        plot(yp, xp, concat([str(colors(rem(kk - 1, 6) + 1)), 'o']))
        plot(calib_data.ocam_model.yc, calib_data.ocam_model.xc, 'ro')
        axis(concat([1, calib_data.width, 1, calib_data.height]))
        drawnow
        set(5 + kk, color=concat([1, 1, 1]))
        set(5 + kk, 'Name', concat(['Image ', str(kk)]), 'NumberTitle', 'off')
        draw_axes(calib_data.Xp_abs[:, :, kk], calib_data.Yp_abs[:, :, kk],
                  calib_data.n_sq_y)
        hold('off')
        zoom('on')
