import os
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import *

import ocam
from ocam import visualize_calibration, reprojectpoints
from ocam.get_checkerboard_corners import get_checkerboard_corners, _refine_corners, \
    _flatten_corner_matrices, _fix_numbering_direction, _fix_nonsquare_board_ambiguity, detect_corners_old

data_dir = Path(os.path.dirname(os.path.realpath(__file__))) / 'data'


@pytest.fixture
def calib_data():
    n_sq_x = 6
    n_sq_y = 5
    dX = 30
    dY = 30
    imgs = sorted(data_dir.glob('*.jpg'))
    return ocam.CalibData(imgs, n_sq_x, n_sq_y, dX, dY)


def reference(title, idx):
    with (data_dir / 'reference' / f'{title}_{idx + 1:02d}.csv').open() as f:
        return np.loadtxt(f, delimiter=',')


@pytest.skip
@pytest.mark.parametrize("kk", np.arange(10))
def test_get_checkerboard_corners(kk, calib_data):
    img_shape = calib_data.height, calib_data.height
    n_sq_x, n_sq_y = calib_data.n_sq_x, calib_data.n_sq_y
    img = calib_data.read_image(kk)

    cornerInfo, cornersX, cornersY = detect_corners_old(calib_data.imgs[kk], n_sq_x, n_sq_y)
    # assert_array_equal(cornersX, reference('x_initial', kk))
    # assert_array_equal(cornersY, reference('y_initial', kk))

    # bounds = _calculate_bounds(cornersX, cornersY, img_shape)
    # assert_array_equal(bounds, reference('bounds', kk))

    cornersX, cornersY = _refine_corners(cornersX, cornersY, img)
    assert_array_almost_equal(np.array(cornersX), reference('x_refined', kk) - 1, decimal=2)
    assert_array_almost_equal(np.array(cornersY), reference('y_refined', kk) - 1, decimal=2)

    x, y = _flatten_corner_matrices(cornersX, cornersY)
    assert_array_almost_equal(np.squeeze(np.array(x)), reference('x_flattened', kk) - 1, decimal=2)
    assert_array_almost_equal(np.squeeze(np.array(y)), reference('y_flattened', kk) - 1, decimal=2)

    x, y = _fix_numbering_direction(x, y, n_sq_x, n_sq_y)
    assert_array_almost_equal(np.squeeze(np.array(x)), reference('x_numbering_direction', kk) - 1, decimal=2)
    assert_array_almost_equal(np.squeeze(np.array(y)), reference('y_numbering_direction', kk) - 1, decimal=2)

    x, y = _fix_nonsquare_board_ambiguity(x, y, n_sq_x, n_sq_y)
    assert_array_almost_equal(np.squeeze(np.array(x)), reference('x_assign_corner', kk) - 1, decimal=2)
    assert_array_almost_equal(np.squeeze(np.array(y)), reference('y_assign_corner', kk) - 1, decimal=2)


def test_get_checkerboard_corners_viz(calib_data):
    get_checkerboard_corners(calib_data, 2, visualize=True)


def test_calibrate(calib_data):
    calib_data.Xp_abs = np.zeros((calib_data.n_imgs, calib_data.n_corners))
    calib_data.Yp_abs = np.zeros((calib_data.n_imgs, calib_data.n_corners))
    for kk in range(10):
        x, y = get_checkerboard_corners(calib_data, kk)
        calib_data.Xp_abs[kk, :] = x
        calib_data.Yp_abs[kk, :] = y
        calib_data.active_images[kk] = True
        calib_data.ima_proc.append(kk)
        calib_data.ima_proc.sort()

    calib_data.taylor_order = 4
    RRfin, ss = ocam.calibrate_linear(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs,
                                      calib_data.Yp_abs, calib_data.ocam_model.xc,
                                      calib_data.ocam_model.yc, calib_data.taylor_order,
                                      calib_data.ima_proc)
    print(ss)
    print(RRfin)
    calib_data.RRfin = RRfin
    calib_data.ocam_model.ss = ss
    reprojectpoints(calib_data)
    visualize_calibration(calib_data)
