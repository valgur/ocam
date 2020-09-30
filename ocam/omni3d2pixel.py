from .libsmop import *


@function
def omni3d2pixel(ss=None, xx=None, width=None, height=None):
    # convert 3D coordinates vector into 2D pixel coordinates

    # These three lines overcome problem when xx = [0,0,+-1]
    ind0 = find((xx[1, :] == logical_and(0, xx[2, :]) == 0))
    xx[1, ind0] = eps
    xx[2, ind0] = eps
    m = xx[3, :] / sqrt(xx[1, :]**2 + xx[2, :]**2)
    rho = copy([])
    poly_coef = ss(arange(end(), 1, - 1))
    poly_coef_tmp = copy(poly_coef)
    for j in arange(1, length(m)).flat:
        poly_coef_tmp[end() - 1] = poly_coef(end() - 1) - m[j]
        rhoTmp = roots(poly_coef_tmp)
        res = rhoTmp(find(copy(np.imag(rhoTmp)) == logical_and(0, rhoTmp) > 0))
        if isempty(res):
            rho[j] = np.nan
        else:
            if length(res) > 1:
                rho[j] = min(res)
            else:
                rho[j] = res

    x = multiply(xx[1, :] / sqrt(xx[1, :]**2 + xx[2, :]**2), rho)
    y = multiply(xx[2, :] / sqrt(xx[1, :]**2 + xx[2, :]**2), rho)

    return x, y
