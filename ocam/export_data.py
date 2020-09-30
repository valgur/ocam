from .calib_data import OCamModel
from .findinvpoly import findinvpoly
from .libsmop import *


@function
def export_data(ocam_model: OCamModel):
    if logical_not(isfield(ocam_model, 'invpol')):
        width = ocam_model.width
        height = ocam_model.height
        ocam_model.invpol = copy(findinvpoly(ocam_model.ss, sqrt((width / 2)**2 + (height / 2)**2)))

    fid = open('calib_results.txt', 'w')
    fprintf(fid,
            '#polynomial coefficients for the DIRECT mapping function (ocam_model.ss in MATLAB). These are used by cam2world\n\n')
    fid.write('%d ' % (length(ocam_model.ss)))

    for i in arange(1, length(ocam_model.ss)).flat:
        fid.write('%e ' % (ocam_model.ss[i]))

    fprintf(fid, '\n\n')
    fprintf(fid,
            '#polynomial coefficients for the inverse mapping function (ocam_model.invpol in MATLAB). These are used by world2cam\n\n')
    fid.write('%d ' % (length(ocam_model.invpol)))

    for i in arange(1, length(ocam_model.invpol)).flat:
        fid.write('%f ' % (ocam_model.invpol(end() - i + 1)))

    fprintf(fid, '\n\n')
    fprintf(fid, '#center: "row" and "column", starting from 0 (C convention)\n\n')
    fid.write('%f %f\n\n' % (ocam_model.xc - 1, ocam_model.yc - 1))
    fprintf(fid, '#affine parameters "c", "d", "e"\n\n')
    fid.write('%f %f %f\n\n' % (ocam_model.c, ocam_model.d, ocam_model.e))
    fprintf(fid, '#image size: "height" and "width"\n\n')
    fid.write('%d %d\n\n' % (ocam_model.height, ocam_model.width))
    fid.close()
