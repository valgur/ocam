# This small script looks in the direcory and checks if the images are there.
# This works only on Matlab 5.x (otherwise, the dir commands works differently)
# (c) Jean-Yves Bouguet - Dec. 27th, 1999

from .calib_data import CalibData
from .libsmop import *


@function
def check_directory(calib_data: CalibData):
    l = dir(concat([calib_data.calib_name, '*']))
    Nl = size(l, 1)
    Nima_valid = 0
    ind_valid = copy([])
    loc_extension = copy([])
    length_name = size(calib_data.calib_name, 2)
    if Nl > 0:
        for pp in arange(1, Nl).flat:
            filenamepp = l(pp).name
            if isempty(calib_data.calib_name):
                iii = 1
            else:
                iii = findstr(filenamepp, calib_data.calib_name)
            loc_ext = findstr(filenamepp, calib_data.format_image)
            string_num = filenamepp(arange(length_name + 1, loc_ext - 2))
            if isempty(float(string_num)):
                iii = copy([])
            if logical_not(isempty(iii)):
                if iii[1] != 1:
                    iii = copy([])
            if logical_and(logical_not(isempty(iii)), logical_not(isempty(loc_ext))):
                Nima_valid += 1
                ind_valid = concat([ind_valid, pp])
                loc_extension = concat([loc_extension, loc_ext[1]])
        if Nima_valid == 0:
            format_image = upper(calib_data.format_image)
            for pp in arange(1, Nl).flat:
                filenamepp = l(pp).name
                if isempty(calib_data.calib_name):
                    iii = 1
                else:
                    iii = findstr(filenamepp, calib_data.calib_name)
                loc_ext = findstr(filenamepp, format_image)
                string_num = filenamepp(arange(length_name + 1, loc_ext - 2))
                if isempty(float(string_num)):
                    iii = copy([])
                if logical_not(isempty(iii)):
                    if iii[1] != 1:
                        iii = copy([])
                if logical_and(logical_not(isempty(iii)), logical_not(isempty(loc_ext))):
                    Nima_valid += 1
                    ind_valid = concat([ind_valid, pp])
                    loc_extension = concat([loc_extension, loc_ext[1]])
            if Nima_valid == 0:
                print('No image found. File format may be wrong.')
            else:
                # Get all the string numbers:
                string_length = zeros(1, Nima_valid)
                indices = zeros(1, Nima_valid)
                for ppp in arange(1, Nima_valid).flat:
                    name = l(ind_valid(ppp)).name
                    string_num = name(arange(length_name + 1, loc_extension(ppp) - 2))
                    string_length[ppp] = size(string_num, 2)
                    indices[ppp] = float(string_num)
                first_num = min(indices)
                n_ima = max(indices) - first_num + 1
                if min(string_length) == max(string_length):
                    N_slots = min(string_length)
                    type_numbering = 1
                else:
                    N_slots = 1
                    type_numbering = 0
                image_numbers = arange(first_num, calib_data.n_imgs - 1 + first_num)
                calib_data.active_images = ones(1, calib_data.n_imgs)
        else:
            # Get all the string numbers:
            string_length = zeros(1, Nima_valid)
            indices = zeros(1, Nima_valid)
            for ppp in arange(1, Nima_valid).flat:
                name = l(ind_valid(ppp)).name
                string_num = name(arange(length_name + 1, loc_extension(ppp) - 2))
                string_length[ppp] = size(string_num, 2)
                indices[ppp] = float(string_num)
            first_num = min(indices)
            calib_data.n_ima = max(indices) - first_num + 1
            if min(string_length) == max(string_length):
                N_slots = min(string_length)
                type_numbering = 1
            else:
                N_slots = 1
                type_numbering = 0
            image_numbers = arange(first_num, calib_data.n_imgs - 1 + first_num)
            calib_data.active_images = ones(1, calib_data.n_imgs)
    else:
        print('No image found. Basename may be wrong.')

    return Nima_valid, image_numbers, type_numbering, N_slots
