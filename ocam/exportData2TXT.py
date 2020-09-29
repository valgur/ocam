from .calib_data import CalibData
from .export_data import export_data
from .libsmop import *


@function
def exportData2TXT(calib_data: CalibData = None):
    if isempty(calib_data.n_imgs) or calib_data.calibrated == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    print('Exporting ocam_model to "calib_results.txt"')
    export_data(calib_data.ocam_model)
    print('done')
