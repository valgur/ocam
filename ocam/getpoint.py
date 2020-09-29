from .libsmop import *


@function
def getpoint(ss, m):
    # Given an image point it returns the 3D coordinates of its correspondent optical ray
    return concat([
        [m[1, :]],
        [m[2, :]],
        [polyval(ss(arange(end(), 1, - 1)), sqrt(m[1, :]**2 + m[2, :]**2))]
    ])
