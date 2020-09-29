from .libsmop import *


@function
def FUNrho(ss, rho):
    return polyval(ss(arange(end(), 1, -1)), rho)
