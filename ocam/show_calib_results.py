from .calib_data import CalibData
from .libsmop import *
from .reprojectpoints_adv import reprojectpoints_adv


@function
def show_calib_results(calib_data: CalibData):
    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    M = concat([calib_data.Xt, calib_data.Yt, zeros(size(calib_data.Xt))])
    reprojectpoints_adv(calib_data.ocam_model, calib_data.RRfin, calib_data.ima_proc, calib_data.Xp_abs,
                        calib_data.Yp_abs, M)
    ss = calib_data.ocam_model.ss
    print(ss)
    xc = calib_data.ocam_model.xc
    print(xc)
    yc = calib_data.ocam_model.yc
    print(yc)
    figure[3]
    set(3, 'Name', 'Calibration results', 'NumberTitle', 'off')
    subplot(2, 1, 1)
    plot(arange(0, floor(calib_data.ocam_model.width / 2)),
         polyval(concat([calib_data.ocam_model.ss(arange(end(), 1, -1))]),
                 concat([arange(0, floor(calib_data.ocam_model.width / 2))])))
    grid('on')
    axis('equal')
    xlabel('Distance \'rho\' from the image center in pixels')
    ylabel('f(rho)')
    title('Forward projection function')

    subplot(2, 1, 2)
    plot(arange(0, floor(calib_data.ocam_model.width / 2)), dot(180 / np.pi, np.arctan2(
        arange(0, floor(calib_data.ocam_model.width / 2)),
        - polyval(concat([calib_data.ocam_model.ss(arange(end(), 1, - 1))]),
                  concat([arange(0, floor(calib_data.ocam_model.width / 2))])))) - 90)
    grid('on')
    xlabel('Distance \'rho\' from the image center in pixels')
    ylabel('Degrees')
    title('Angle of optical ray as a function of distance from circle center (pixels)')
