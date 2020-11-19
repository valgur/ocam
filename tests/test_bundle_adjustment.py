import os
from pathlib import Path

import pytest

import ocam
import numpy as np

data_dir = Path(os.path.dirname(os.path.realpath(__file__))) / 'data'


@pytest.fixture
def calib_data():
    n_sq_x = 6
    n_sq_y = 5
    dX = 30
    dY = 30
    imgs = sorted(data_dir.glob('*.jpg'))
    return ocam.CalibData(imgs, n_sq_x, n_sq_y, dX, dY)


def test_bundle_err(calib_data):
    ocam.extract_corners(calib_data, visualize=False)
    ocam.calibration(calib_data, visualize=False)
    ocam.bundle_adjustment(calib_data, robust=True)
    print(calib_data.ocam_model)
    print(calib_data.RRfin)
    r = np.sqrt(calib_data.ocam_model.width**2 + calib_data.ocam_model.height**2) / 2
    inv_poly_res = ocam.findinvpoly(calib_data.ocam_model.ss, r)
    print(inv_poly_res)
