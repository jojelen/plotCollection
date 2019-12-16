"""
Functions for x,y,z data with bins in x,y
"""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, MultipleLocator, FormatStrFormatter
import matplotlib

import numpy as np
import os

###############################################################################
# Utils
###############################################################################


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

###############################################################################
# Binning matrices that pick out max/min/mean of z values
###############################################################################


def meanMatrix(x, y, z, bins, nPad, xmax, xmin, ymax, ymin, zmax, zmin):
    """Returns a bin matrix

    Given the arrays x,y,z, it creates a matrix  with the mean z value in each 
    bin."""

    xWidth = (xmax - xmin) / bins
    yWidth = (ymax - ymin) / bins

    nMatrix = np.zeros((bins, bins))
    totMatrix = np.zeros((bins, bins))

    # (nPad sets number of zero bins in edge)
    averageMatrix = np.zeros((bins + 2 * nPad, bins + 2 * nPad))

    # First, calculate number of points in each bin.
    # We also sum up their z value
    for i in range(len(x)):
        if (x[i] > xmax + nPad*xWidth or x[i] < xmin - nPad*xWidth):
            continue
        if (y[i] > ymax + nPad*yWidth or y[i] < ymin - nPad*yWidth):
            continue
        if (z[i] > zmax or z[i] < zmin):
            continue
        try:
            xbin = int((x[i] - xmin) / xWidth)
            if (xbin == bins):
                xbin = bins - 1
        except:
            print('[ERROR]: For x[i] = {}, xmin = {}, xWidth = {}'.format(x[i],
                                                                          xmin, xWidth))
            raise
        ybin = int((y[i] - ymin) / yWidth)
        if (ybin == bins):
            ybin = bins - 1
        # The imshow() Plots the matrix as an image.
        # Therefore, we reverse the y-axis and transposes the
        # matrix.
        nMatrix[bins - 1 - ybin, xbin] += 1
        totMatrix[bins - 1 - ybin, xbin] += z[i]

    # Calculate average z value
    for i in range(bins):
        for j in range(bins):
            if (nMatrix[i, j] != 0):
                averageMatrix[i + nPad, j +
                              nPad] = totMatrix[i, j] / nMatrix[i, j]

    return averageMatrix


def maxMatrix(x, y, z, bins, nPad, xmax, xmin, ymax, ymin, zmax, zmin):
    """Returns a bin matrix

    Given the arrays x,y,z, it creates a matrix  with the max z value in each 
    bin."""

    xWidth = (xmax - xmin) / bins
    yWidth = (ymax - ymin) / bins

    # (nPad sets number of zero bins in edge)
    finalMatrix = np.zeros((bins + 2 * nPad, bins + 2 * nPad))

    # Go through each parameter point.
    # Calculate thecorresponding bin.
    # Set bin value to z value if it is larger.
    for i in range(len(x)):
        if (x[i] > xmax + nPad*xWidth or x[i] < xmin - nPad*xWidth):
            continue
        if (y[i] > ymax + nPad*yWidth or y[i] < ymin - nPad*yWidth):
            continue
        if (z[i] > zmax or z[i] < zmin):
            continue
        try:
            xbin = int((x[i] - xmin) / xWidth)
            if (xbin == bins):
                xbin = bins - 1
        except:
            print('[ERROR]: For x[i] = {}, xmin = {}, xWidth = {}, xmin = {}, xmax = {}'.format(x[i],
                                                                                                xmin, xWidth, xmin, xmax))
            raise
        ybin = int((y[i] - ymin) / yWidth)
        if (ybin == bins):
            ybin = bins - 1
        # The imshow() Plots the matrix as an image.
        # Therefore, we reverse the y-axis and transposes the
        # matrix.
        if (z[i] > finalMatrix[bins - 1 - ybin + nPad, xbin + nPad]):
            finalMatrix[bins - 1 - ybin + nPad, xbin + nPad] = z[i]

    return finalMatrix


def minMatrix(x, y, z, bins, nPad, xmax, xmin, ymax, ymin, zmax, zmin):
    """Returns a bin matrix

    Given the arrays x,y,z, it creates a matrix  with the min z value in each 
    bin."""

    xWidth = (xmax - xmin) / bins
    yWidth = (ymax - ymin) / bins

    # (nPad sets number of zero bins in edge)
    # A matrix with the maximum z values
    finalMatrix = maxMatrix(x, y, z, bins, nPad, xmax,
                            xmin, ymax, ymin, zmax, zmin)

    # Go through each parameter point.
    # Calculate thecorresponding bin.
    # Set bin value to z value if it is smaller.
    for i in range(len(x)):
        if (x[i] > xmax + nPad*xWidth or x[i] < xmin - nPad*xWidth):
            continue
        if (y[i] > ymax + nPad*yWidth or y[i] < ymin - nPad*yWidth):
            continue
        if (z[i] > zmax or z[i] < zmin):
            continue
        xbin = int((x[i] - xmin) / xWidth)
        if (xbin == bins):
            xbin = bins - 1
        ybin = int((y[i] - ymin) / yWidth)
        if (ybin == bins):
            ybin = bins - 1
        # The imshow() Plots the matrix as an image.
        # Therefore, we reverse the y-axis and transposes the
        # matrix.
        if (z[i] < finalMatrix[bins - 1 - ybin + nPad, xbin + nPad] and z[i] > 0):
            finalMatrix[bins - 1 - ybin + nPad, xbin + nPad] = z[i]

    return finalMatrix

###############################################################################
# Plot function
###############################################################################


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
        norming=LogNorm(),
        output=None):
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
        xyMatrix = meanMatrix(x, y, z, bins, nPad, xmax,
                              xmin, ymax, ymin, zmax, zmin)
    elif (measure == 'max'):
        xyMatrix = maxMatrix(x, y, z, bins, nPad, xmax,
                             xmin, ymax, ymin, zmax, zmin)
    elif (measure == 'min'):
        xyMatrix = minMatrix(x, y, z, bins, nPad, xmax,
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
    if output == None:
        plt.show()
    else:
        try:
            directory = os.path.dirname(output)
            if not os.path.exists(directory):
                os.makedirs(directory)
            plt.savefig(output)
        except Exception as e:
            print('[ERROR]: Could not save {}\n\t {}'.format(output, e))
