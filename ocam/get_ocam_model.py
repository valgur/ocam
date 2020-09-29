from .libsmop import *


@function
def get_ocam_model(radius):
    ss = concat(
        [[- 140.5116937602191], [0], [0.0002716608082380784], [5.257341861497706e-06], [- 1.067888507955045e-09]])
    k = radius / 498.0
    kk = dot(ones(size(ss)), k)
    kkk = kk**(concat([arange(- 1, length(ss) - 2)]).T)
    ss_n = ss / kkk
    # theta_cam=rad2deg(np.arctan2(polyval([ss_n[end:-1:1]],[0:floor(ny/2)]),[0:floor(ny/2)]));
    return ss_n
