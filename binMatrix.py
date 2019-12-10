"""
Functions for creating 2D grids with bins
"""
import numpy as np

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
    finalMatrix = maxMatrix(x, y, z, bins, nPad, xmax,xmin,ymax,ymin,zmax,zmin)

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
