from .libsmop import *


@function
def planefrompoints(xp, yp, zp):
    A = concat([xp, yp, zp])
    B = -ones(size(xp))
    return pinv(A) @ B
