from .calib_data import CalibData
from .libsmop import *


@function
def saving_calib(calib_data: CalibData):
    if logical_or(isempty(calib_data.n_imgs), calib_data.calibrated) == 0:
        print(
            '\nNo calibration data available. You must first calibrate your camera.\nClick on "Calibration" or "Find center"\n')
        return

    save_name = 'Omni_Calib_Results'
    if exist(concat([save_name, '.mat'])) == 2:
        disp('WARNING: File Omni_Calib_Results.mat already exists')
        if exist('copyfile'):
            pfn = - 1
            cont = 1
            while cont:
                pfn += 1
                postfix = concat(['_old', str(pfn)])
                save_name = concat(['Omni_Calib_Results', postfix])
                cont = (exist(concat([save_name, '.mat'])) == 2)

            copyfile('Omni_Calib_Results.mat', concat([save_name, '.mat']))
            disp(concat(['Copying the current Omni_Calib_Results.mat file to ', save_name, '.mat']))
            if exist('Omni_Calib_Results.m') == 2:
                copyfile('Omni_Calib_Results.m', concat([save_name, '.m']))
                disp(concat(['Copying the current Omni_Calib_Results.m file to ', save_name, '.m']))
            cont_save = 1
        else:
            disp('The file Omni_Calib_Results.mat is about to be changed.')
            cont_save = input('Do you want to continue? ([]=no,other=yes) ', 's')
            cont_save = logical_not(isempty(cont_save))
    else:
        cont_save = 1

    if cont_save:
        if exist('calib_name'):
            print(concat(['\nSaving calibration results under ', save_name, '.mat, please wait ...\n']), end='')
        save_name = 'Omni_Calib_Results'
        data_name = 'calib_data'
        string_save = concat(['save ', save_name, ' ', data_name])
        eval(string_save)
        print('done')
