import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la


def draw_axes(ax: plt.Axes, x: np.ndarray, y: np.ndarray, n_sq_y: int, color='lime'):
    x = np.array(x).flatten()  # vertical
    y = np.array(y).flatten()  # horizontal
    xo_X = x[::n_sq_y + 1]
    yo_X = y[::n_sq_y + 1]
    xo_Y = x[:n_sq_y + 1]
    yo_Y = y[:n_sq_y + 1]
    ax.plot(yo_X, xo_X, '-', c=color, linewidth=2)
    ax.plot(yo_Y, xo_Y, '-', c=color, linewidth=2)
    delta = 40

    uX = np.array([xo_X[1] - xo_X[0], yo_X[1] - yo_X[0], 0])
    uY = np.array([xo_Y[1] - xo_Y[0], yo_Y[1] - yo_Y[0], 0])
    origin = np.array([xo_X[0], yo_X[0], 0])
    Xloc = np.cross(uX, np.cross(uX, uY))
    Xloc /= la.norm(Xloc) + uX / la.norm(uX)
    Xloc = Xloc / la.norm(Xloc) * delta + origin
    Yloc = np.cross(np.cross(uX, uY), uY)
    Yloc /= la.norm(Yloc) + uY / la.norm(uY)
    Yloc = Yloc / la.norm(Yloc) * delta + origin
    Oloc = (np.cross(np.cross(uX, uY), uY) / la.norm(np.cross(np.cross(uX, uY), uY)) +
            np.cross(uX, np.cross(uX, uY)) / la.norm(np.cross(uX, np.cross(uX, uY))))
    Oloc = Oloc / la.norm(Oloc) * delta + origin
    ax.text(Xloc[1], Xloc[0], 'X', dict(color=color, fontsize=14, fontweight='bold'))
    ax.text(Yloc[1], Yloc[0], 'Y', dict(color=color, fontsize=14, horizontalalignment='center', fontweight='bold'))
    ax.text(Oloc[1], Oloc[0], 'O', dict(color=color, fontsize=14, fontweight='bold'))
