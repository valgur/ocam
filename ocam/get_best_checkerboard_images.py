from .get_checkerboard_corners import get_checkerboard_corners
from .libsmop import *


@function
def get_best_checkerboard_images(ima_numbers=None, I=None, n_sq_x=None, n_sq_y=None, use_corner_find=None):
    count = 0
    for kk in ima_numbers.flat:
        callBack, x, y = get_checkerboard_corners(I, kk, n_sq_x, n_sq_y, use_corner_find, nargout=3)
        if callBack > 0:
            count += 1
            ima_proc[count] = kk
            Xp_abs[:, :, kk] = x
            Yp_abs[:, :, kk] = y
    return ima_proc, Xp_abs, Yp_abs
