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
import os
from dataclasses import InitVar, dataclass, field

import numpy as np
from imageio import imread


@dataclass
class OCamModel:
    img_shape: InitVar[(int, int)] = None  # (height, width) for initialization

    ss: np.ndarray = None  # Coefficients of polynomial
    xc: float = None  # x-coordinate of image center
    yc: float = None  # y-coordinate of image center
    # Affine transformation parameters
    c: float = 1
    d: float = 0
    e: float = 0
    width: int = None  # Image width
    height: int = None  # Image height
    pol: np.ndarray = None  # Inverse polynomial (optional)

    def __post_init__(self, img_shape):
        if img_shape:
            self.height = img_shape[0]
            self.width = img_shape[1]
            self.xc = round(self.height / 2.)
            self.yc = round(self.width / 2.)


@dataclass()
class CalibData:
    """Stores image and calibration data used by the calibration toolbox"""

    imgs: list  # Image file names
    n_sq_x: int  # Number of squares in x-direction
    n_sq_y: int  # Number of squares in y-direction
    dX: float  # Width of a square on checkerboard (mm)
    dY: float  # Height of a square on checkerboard (mm)

    taylor_order: int = 4  # Order of the polynomial

    active_images: np.ndarray = None  # Vector indicating images used
    ima_proc: list = field(default_factory=list)  # Images being processed

    # Calibration model
    ocam_model: OCamModel = OCamModel()

    RRfin: np.ndarray = None  # Extrinsic parameters of checkerboards

    Xt: np.ndarray = None  # Checkerboard corner coordinates (mm)
    Yt: np.ndarray = None
    Xp_abs: np.ndarray = None  # Checkerboard corner coordinates (px)
    Yp_abs: np.ndarray = None

    wintx: int = None  # Size of corner search window for assisted
    winty: int = None  # manual corner selection

    # Flags
    no_image_file: bool = None  # Indicates missing image files
    calibrated: bool = False  # Indicates calibrated ocam_model
    optimized: bool = False  # Non-linear optimization has been applied

    map: str = 'gray'  # Colormap for displaying images

    # Non-linear optimization result from scipy.optimize.least_squares
    ba_result: dict = None
    # weights from robust optimization
    weights: np.ndarray = None

    @property
    def n_imgs(self):
        """Number of input images"""
        return len(self.imgs)

    @property
    def ind_active(self):
        """Indices of images used"""
        return np.argwhere(self.active_images).flatten()

    @property
    def width(self):
        return self.ocam_model.width

    @property
    def height(self):
        return self.ocam_model.height

    @property
    def n_cor_x(self):
        return self.n_sq_x + 1

    @property
    def n_cor_y(self):
        return self.n_sq_y + 1

    @property
    def n_corners(self):
        return self.n_cor_x * self.n_cor_y

    def __post_init__(self):
        if self.n_imgs == 0:
            raise ValueError('No input images provided')

        for f in self.imgs:
            if not os.path.exists(f) and os.path.isfile(f):
                raise ValueError(f"Image '{f}' does not exist")

        if self.active_images is None or len(self.active_images) == 0:
            self.active_images = np.ones(self.n_imgs, dtype=bool)
        else:
            assert len(self.active_images) == self.n_imgs

        img = self.read_image(0)
        self.ocam_model = OCamModel(img.shape[:2])

        # Arranging the pixel of the world
        Xt = []
        Yt = []
        for i in range(self.n_sq_x + 1):
            for j in range(self.n_sq_y + 1):
                Yt.append(j * self.dY)
                Xt.append(i * self.dX)
        self.Xt = np.array(Xt)
        self.Yt = np.array(Yt)

    def read_image(self, idx):
        img = imread(self.imgs[idx])
        if len(img.shape) > 2:
            img = img @ np.array([0.299, 0.587, 0.114])
        return img
