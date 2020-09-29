from .calib_data import CalibData
from .libsmop import *


@function
def mosaic(calib_data: CalibData):
    if isempty(calib_data.I):
        active_images_save = calib_data.active_images
        ima_read_calib(calib_data)
        calib_data.active_images = copy(active_images_save)

    if len(calib_data.ind_read) == 0:
        return

    n_col = floor(sqrt(calib_data.n_imgs * calib_data.ocam_model.width / calib_data.ocam_model.height))
    n_row = ceil(calib_data.n_imgs / n_col)
    ker2 = 1
    for ii in arange(1, n_col).flat:
        ker2 = conv(ker2, concat([1 / 4, 1 / 2, 1 / 4]))

    II = calib_data.I[1](arange(1, end(), n_col), arange(1, end(), n_col))
    ny2, nx2 = size(II, nargout=2)
    kk_c = 1
    II_mosaic = copy([])
    for jj in arange(1, n_row).flat:
        II_row = copy([])
        for ii in arange(1, n_col).flat:
            if logical_and((kk_c <= calib_data.n_imgs), logical_not(isempty(calib_data.I[kk_c]))):
                if calib_data.active_images(kk_c):
                    I = calib_data.I[kk_c]
                    I = I(arange(1, end(), n_col), arange(1, end(), n_col))
                else:
                    I = zeros(ny2, nx2)
            else:
                I = zeros(ny2, nx2)
            II_row = concat([II_row, I])
            if ii != n_col:
                II_row = concat([II_row, zeros(ny2, 3)])
            kk_c += 1
        nn2 = size(II_row, 2)
        if jj != n_row:
            II_row = concat([[II_row], [zeros(3, nn2)]])
        II_mosaic = concat([[II_mosaic], [II_row]])

    figure[2]
    image(II_mosaic)
    colormap(gray[256])
    title('Calibration images')
    set(gca, 'Xtick', [])
    set(gca, 'Ytick', [])
    axis('image')
