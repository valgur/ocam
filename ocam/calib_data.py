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
# filename: C_self.m
# Code was added to 
# save statistical relevant variables

from dataclasses import dataclass

from imageio import imread

from .libsmop import *


@dataclass
class OCamModel:
    ss: matlabarray = None  # Coefficients of polynomial
    xc: float = None  # x-coordinate of image center
    yc: float = None  # y-coordinate of image center
    # Affine transformation parameters
    c: float = None
    d: float = None
    e: float = None
    width: int = None  # Image width
    height: int = None  # Image height


@dataclass
class StatEO:
    stdEO: matlabarray = None
    varEO: matlabarray = None
    sg0: matlabarray = None
    Exx: matlabarray = None


@dataclass
class StatIO:
    stdIO: matlabarray = None
    varIO: matlabarray = None
    sg0: matlabarray = None
    Exx: matlabarray = None


@dataclass()
class CalibData:
    """Stores image and calibration data used by the calibration toolbox"""

    imgs: matlabarray  # Image file names
    n_sq_x: int  # Number of squares in x-direction
    n_sq_y: int  # Number of squares in y-direction
    dX: float  # Width of a square on checkerboard (mm)
    dY: float  # Height of a square on checkerboard (mm)

    taylor_order: int = 4  # Order of the polynomial

    active_images: matlabarray = None  # Vector indicating images used
    ima_proc: matlabarray = None  # Images being processed

    # Calibration model
    ocam_model: OCamModel = OCamModel()

    RRfin: matlabarray = None  # Extrinsic parameters of checkerboards

    Xt: float = None  # Checkerboard corner coordinates (mm)
    Yt: float = None
    Xp_abs: float = None  # Checkerboard corner coordinates (px)
    Yp_abs: float = None
    Xp_abss: float = None  # Checkerboard sub pixel corner coordinates (px)
    Yp_abss: float = None

    wintx: int = None  # Size of corner search window for assisted
    winty: int = None  # manual corner selection

    # Flags
    no_image_file: bool = None  # Indicates missing image files
    calibrated: bool = False  # Indicates calibrated ocam_model

    map: str = 'gray'  # Colormap for displaying images

    ## ================
    #  added code
    #
    # errors
    # overall errors
    errMean: float = None
    errStd: float = None
    mse: float = None
    rms: float = None
    runtime: float = None
    rmsAfterCenter: float = None

    statEO: StatEO = StatEO()
    statIO: StatIO = StatIO()

    # weights from robust optimization
    weights: matlabarray = None

    optimized: bool = False

    ## ================

    @property
    def n_imgs(self):
        """Number of input images"""
        return len(self.imgs)

    @property
    def ind_active(self):
        """Indices of images used"""
        return find(self.active_images)

    def __post_init__(self):
        self.imgs = copy(self.imgs)
        if self.n_imgs == 0:
            raise ValueError('No input images provided')

        for f in self.imgs:
            if not os.path.exists(f) and os.path.isfile(f):
                raise ValueError(f"Image '{f}' does not exist")

        if self.active_images is None or len(self.active_images) == 0:
            self.active_images = ones(1, self.n_imgs)
        else:
            assert len(self.active_images) == self.n_imgs

        self.Xt = copy([])
        self.Yt = copy([])
        self.ocam_model.xc = round(self.ocam_model.height / 2.)
        self.ocam_model.yc = round(self.ocam_model.width / 2.)

        img = self.read_image(1)
        self.ocam_model.height, self.ocam_model.width = img.shape[:2]

    def read_image(self, idx):
        img = imread(self.imgs[idx])
        if img.shape[2] > 1:
            img = 0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2]
        img = matlabarray(np.squeeze(img))
        return img
