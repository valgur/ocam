import os
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import *

import ocam
from ocam import call_find_corners
from ocam.get_checkerboard_corners import get_checkerboard_corners, _pixel_coords_to_matlab_format, _calculate_bounds, \
    _refine_corners, _flatten_corner_matrices, _fix_numbering_direction, _fix_nonsquare_board_ambiguity

@pytest.fixture
def data_dir():
    return Path(os.path.dirname(os.path.realpath(__file__))) / 'data'


@pytest.fixture
def default_calib_data(data_dir):
    n_sq_x = 6
    n_sq_y = 5
    dX = 30
    dY = 30
    imgs = sorted(data_dir.glob('*.jpg'))
    return ocam.CalibData(imgs, n_sq_x, n_sq_y, dX, dY)


expected_cornersX_in = np.asarray([
    [677, 705.5, 739, 776, 816, 858.5],
    [671, 701, 735.5, 773.5, 814.5, 859],
    [661, 691, 725.5, 764.5, 807, 852],
    [648, 677, 711, 749, 791, 834.5],
    [632, 658.5, 690, 725.5, 764.5, 806],
    [614, 637.5, 665, 696.5, 730.5, 768],
    [595.5, 614.5, 638, 663.5, 693, 725]
])
expected_cornersY_in = np.asarray([
    [372.5, 374, 378, 383, 388.5, 397],
    [395, 400, 409, 417.5, 429.5, 442],
    [419.5, 430.5, 442.5, 457, 473, 491],
    [446, 460, 478, 497.5, 520, 543.5],
    [471.5, 490, 512, 536.5, 564.5, 592.5],
    [494.5, 516.5, 542.5, 571, 602, 634.5],
    [516, 540, 567, 597.5, 631.5, 666.5]
])


def test_get_checkerboard_corners(default_calib_data):
    calib_data = default_calib_data
    kk = 1
    img_shape = calib_data.height, calib_data.height
    n_sq_x, n_sq_y = calib_data.n_sq_x, calib_data.n_sq_y
    img = calib_data.read_image(kk)

    cornerInfo, cornersX, cornersY = call_find_corners(calib_data.imgs[kk - 1], n_sq_x, n_sq_y)
    assert_array_equal(cornersX, expected_cornersX_in)
    assert_array_equal(cornersY, expected_cornersY_in)

    cornersX, cornersY = _pixel_coords_to_matlab_format(cornersX, cornersY)
    assert_array_equal(cornersX, expected_cornersX_in + 1)
    assert_array_equal(cornersY, expected_cornersY_in + 1)

    bounds = _calculate_bounds(cornersX, cornersY, img_shape)
    assert_array_equal(bounds, [543.8, 768, 314.7, 738.06])

    cornersX, cornersY = _refine_corners(cornersX, cornersY, img)
    assert_array_almost_equal(np.array(cornersX), np.array([
        [678.69, 707.52, 740.26, 777.03, 816.78, 859.2],
        [671.78, 701.23, 736.1, 774.3, 816.16, 860.14],
        [662.2, 692.02, 726.74, 765.61, 808.41, 853.1],
        [649.18, 678.2, 711.75, 749.71, 791.74, 835.68],
        [633.41, 659.77, 691.02, 726.41, 765.78, 807.33],
        [615.09, 638.48, 665.87, 696.81, 731.87, 768.83],
        [596.54, 615.94, 638.58, 664.77, 694.34, 726.36]
    ]), decimal=2)
    assert_array_almost_equal(np.array(cornersY), np.array([
        [372.75, 375.42, 379.22, 384.15, 390.25, 397.5],
        [394.98, 401.44, 409.49, 419.2, 430.4, 442.94],
        [420.6, 431.16, 443.54, 457.93, 474.42, 492.24],
        [446.66, 461.51, 478.76, 498.61, 520.73, 544.23],
        [472.11, 490.94, 513.06, 537.69, 565.42, 594.11],
        [495.7, 517.92, 543.51, 571.81, 603.2, 635.6],
        [517.12, 540.92, 568.44, 598.83, 632.55, 667.41]
    ]), decimal=2)

    x, y = _flatten_corner_matrices(cornersX, cornersY)
    assert_array_almost_equal(
        np.array(x),
        [[678.69, 707.52, 740.26, 777.03, 816.78, 859.2, 671.78, 701.23, 736.1, 774.3, 816.16, 860.14, 662.2, 692.02,
          726.74, 765.61, 808.41, 853.1, 649.18, 678.2, 711.75, 749.71, 791.74, 835.68, 633.41, 659.77, 691.02,
          726.41, 765.78, 807.33, 615.09, 638.48, 665.87, 696.81, 731.87, 768.83, 596.54, 615.94, 638.58, 664.77,
          694.34, 726.36]], decimal=2)
    assert_array_almost_equal(
        np.array(y),
        [[372.75, 375.42, 379.22, 384.15, 390.25, 397.5, 394.98, 401.44, 409.49, 419.2, 430.4, 442.94, 420.6, 431.16,
          443.54, 457.93, 474.42, 492.24, 446.66, 461.51, 478.76, 498.61, 520.73, 544.23, 472.11, 490.94, 513.06,
          537.69, 565.42, 594.11, 495.7, 517.92, 543.51, 571.81, 603.2, 635.6, 517.12, 540.92, 568.44, 598.83, 632.55,
          667.41]], decimal=2)

    x, y = _fix_numbering_direction(x, y, n_sq_x, n_sq_y)
    assert_array_almost_equal(
        np.array(x),
        [[678.69, 707.52, 740.26, 777.03, 816.78, 859.2, 671.78, 701.23, 736.1, 774.3, 816.16, 860.14, 662.2, 692.02,
          726.74, 765.61, 808.41, 853.1, 649.18, 678.2, 711.75, 749.71, 791.74, 835.68, 633.41, 659.77, 691.02, 726.41,
          765.78, 807.33, 615.09, 638.48, 665.87, 696.81, 731.87, 768.83, 596.54, 615.94, 638.58, 664.77, 694.34,
          726.36]], decimal=2)
    assert_array_almost_equal(
        np.array(y),
        [[372.75, 375.42, 379.22, 384.15, 390.25, 397.5, 394.98, 401.44, 409.49, 419.2, 430.4, 442.94, 420.6, 431.16,
          443.54, 457.93, 474.42, 492.24, 446.66, 461.51, 478.76, 498.61, 520.73, 544.23, 472.11, 490.94, 513.06,
          537.69,
          565.42, 594.11, 495.7, 517.92, 543.51, 571.81, 603.2, 635.6, 517.12, 540.92, 568.44, 598.83, 632.55, 667.41]],
        decimal=2)

    x, y = _fix_nonsquare_board_ambiguity(x, y, n_sq_x, n_sq_y)
    assert_array_almost_equal(
        np.array(x),
        [[678.69, 707.52, 740.26, 777.03, 816.78, 859.2, 671.78, 701.23, 736.1, 774.3, 816.16, 860.14, 662.2, 692.02,
          726.74, 765.61, 808.41, 853.1, 649.18, 678.2, 711.75, 749.71, 791.74, 835.68, 633.41, 659.77, 691.02, 726.41,
          765.78, 807.33, 615.09, 638.48, 665.87, 696.81, 731.87, 768.83, 596.54, 615.94, 638.58, 664.77, 694.34,
          726.36]], decimal=2)
    assert_array_almost_equal(
        np.array(y),
        [[372.75, 375.42, 379.22, 384.15, 390.25, 397.5, 394.98, 401.44, 409.49, 419.2, 430.4, 442.94, 420.6, 431.16,
          443.54, 457.93, 474.42, 492.24, 446.66, 461.51, 478.76, 498.61, 520.73, 544.23, 472.11, 490.94, 513.06,
          537.69,
          565.42, 594.11, 495.7, 517.92, 543.51, 571.81, 603.2, 635.6, 517.12, 540.92, 568.44, 598.83, 632.55, 667.41]],
        decimal=2)


def test_get_checkerboard_corners_viz(default_calib_data):
    get_checkerboard_corners(default_calib_data, 1, refine_corners=True, visualize=True)
