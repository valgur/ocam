import matplotlib.pyplot as plt
import numpy as np

from .calib_data import CalibData
from .calibrate import calibrate_linear
from .get_checkerboard_corners import get_checkerboard_corners
from .reprojectpoints import reprojectpoints


def extract_corners(calib_data: CalibData, visualize=True):
    if calib_data.Xp_abs is None:
        calib_data.Xp_abs = np.zeros((calib_data.n_imgs, calib_data.n_corners))
        calib_data.Yp_abs = np.zeros((calib_data.n_imgs, calib_data.n_corners))
    count = 0
    for kk in calib_data.ind_active:
        if kk in calib_data.ima_proc:
            continue
        x, y = get_checkerboard_corners(calib_data, kk, visualize=visualize)
        if x is None:
            # Automatic corner extraction failed
            continue
        count += 1
        calib_data.Xp_abs[kk, :] = x
        calib_data.Yp_abs[kk, :] = y
        calib_data.active_images[kk] = True
        calib_data.ima_proc.append(kk)
        calib_data.ima_proc.sort()


def calibration(calib_data: CalibData, visualize=True):
    if not calib_data.ima_proc:
        raise ValueError('No corner data available. Extract grid corners before calibrating.')

    calib_data.calibrated = False
    calib_data.ocam_model.c = 1
    calib_data.ocam_model.d = 0
    calib_data.ocam_model.e = 0
    calib_data.RRfin, calib_data.ocam_model.ss = calibrate_linear(calib_data.Xt, calib_data.Yt, calib_data.Xp_abs,
                                                                  calib_data.Yp_abs, calib_data.ocam_model.xc,
                                                                  calib_data.ocam_model.yc, calib_data.taylor_order,
                                                                  calib_data.ima_proc)
    calib_data.calibrated = True

    print(calib_data.ocam_model.ss)

    reprojectpoints(calib_data)
    if visualize:
        visualize_calibration(calib_data)


def visualize_calibration(calib_data):
    ss = calib_data.ocam_model.ss
    fig, ax = plt.subplots(2)
    fig.suptitle('Calibration results')
    rho = np.arange(int(np.linalg.norm([calib_data.width, calib_data.height])) // 2)
    ax[0].plot(rho, np.polyval(ss[::-1], rho))
    ax[0].grid()
    # ax[0].set_aspect(1)
    ax[0].set_xlabel('Distance \'rho\' from the image center in pixels')
    ax[0].set_ylabel('f(rho)')
    ax[0].set_title('Forward projection function')

    ax[1].plot(rho, 180 / np.pi * np.arctan2(rho, -np.polyval(ss[::-1], rho)) - 90)
    ax[1].grid()
    ax[1].set_xlabel('Distance \'rho\' from the image center in pixels')
    ax[1].set_ylabel('Degrees')
    ax[1].set_title('Angle of optical ray as a function of distance from circle center (pixels)')
    fig.show()
