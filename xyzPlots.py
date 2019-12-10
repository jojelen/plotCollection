import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, MultipleLocator, FormatStrFormatter
import matplotlib

import numpy as np
import binMatrix


def setPlotStyle():
    """Sets font of labels and ticks"""
    # GENERAL SETTINGS
    # For more setting see matplotlib.rcParams.keys()

    # The axes attributes need to be set before the call to subplot
    # Use , weight='bold' for boldface
    plt.rc('font', weight='normal', size=26)
    plt.rc('ytick.major', size=10, pad=6)
    plt.rc('ytick.minor', size=6, pad=6)
    plt.rc('xtick.major', size=10, pad=6)
    plt.rc('xtick.minor', size=6, pad=6)
    plt.rc('xtick', labelsize=26)
    plt.rc('ytick', labelsize=26)
    plt.rc('xtick', top=True)
    plt.rc('ytick', right=True)

    # Customize matplotlib
    matplotlib.rcParams.update(
        {
            'font.sans-serif': 'cmr10',  # the default sans-serif font is cmss10
            # 'font.family': 'sans-serif', #  has trouble with minus signs
            'font.family': 'stixgeneral',
            'mathtext.fontset': 'cm',
            # 'mathtext.fontset': 'dejavuserif', #little bit more 'fancy' style
            'text.usetex': False,
            'axes.unicode_minus': False,  # bug fix related to minus signs
            'axes.linewidth': 1.5
        }
    )
    return


def Heatmap(
        x,
        y,
        z,
        xLabel='x',
        yLabel='y',
        zLabel='z',
        measure='max',
        bins=10,
        nPad=1,
        xmax=None,
        xmin=None,
        ymax=None,
        ymin=None,
        zmax=None,
        zmin=None,
        norming=LogNorm()):
    """Heatmap of three vectors x, y, z"""

    assert len(x) == len(y), "x and y lengths do not match"
    assert len(x) == len(z), "x and z lengths do not match"

    if xmin == None:
        xmin = x.min()
    if xmax == None:
        xmax = x.max()
    if ymin == None:
        ymin = y.min()
    if ymax == None:
        ymax = y.max()
    if zmin == None:
        zmin = z.min()
    if zmax == None:
        zmax = z.max()
    # Width of bins
    xWidth = (xmax - xmin) / bins
    yWidth = (ymax - ymin) / bins

    # Creates a matrix of all bins with a z value in each one
    if (measure == 'mean'):
        xyMatrix = binMatrix.meanMatrix(x, y, z, bins, nPad, xmax,
                                        xmin, ymax, ymin, zmax, zmin)
    elif (measure == 'max'):
        xyMatrix = binMatrix.maxMatrix(x, y, z, bins, nPad, xmax,
                                       xmin, ymax, ymin, zmax, zmin)
    elif (measure == 'min'):
        xyMatrix = binMatrix.minMatrix(x, y, z, bins, nPad, xmax,
                                       xmin, ymax, ymin, zmax, zmin)
    else:
        print('Invalid measure in xyzHeatmap!')
        return

    # Start plotting
    plt.clf()
    fig = plt.figure(figsize=(7, 6))
    setPlotStyle()

    # Colormaps the empty bins as white
    cMap = matplotlib.cm.gnuplot
    cMap.set_bad(color='white')  # Color where heatmap = 0
    xyMatrix = np.ma.masked_where(xyMatrix == 0, xyMatrix)

    plt.xlabel(xLabel)
    plt.ylabel(yLabel)

    plt.grid(alpha=0.5)

    # Plotting the matrix with a z value in each bin
    plt.imshow(xyMatrix,
               cmap=cMap,
               norm=norming,
               extent=[xmin - xWidth * nPad, xmax + xWidth * nPad,
                       ymin - yWidth * nPad, ymax + yWidth * nPad],
               vmin=zmin, vmax=zmax, aspect='auto')

    # Adding a colorbar
    cbar = plt.colorbar(pad=0.05)
    cbar.set_clim(zmin, zmax)
    cbar.set_label(zLabel)
    # Removes minor ticks if there are too many major ticks
    cbar.set_ticks(LogLocator())

    # Ticks configuration
    ax = plt.gca()  # gca returns current axis
    ax.tick_params(
        axis='both',
        which='both',
        direction='in')  # Setting ticks inwards

    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    xMajorTickSpace = plt.xticks()[(0)][1] - plt.xticks()[(0)][0]
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xMajorTickSpace / 4))
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    yMajorTickSpace = plt.yticks()[(0)][1] - plt.yticks()[(0)][0]
    ax.yaxis.set_minor_locator(plt.MultipleLocator(yMajorTickSpace / 4))

    plt.tight_layout()
    plt.show()
