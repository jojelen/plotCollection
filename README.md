# Plot collection

Creates heatmap plots from 3 dimensional data organized in 3 vectors, x,y and z.

## Requirements

matplotlib, numpy

## Example

x,y,z can be a collection of parameter points (z being the temperature across
x,y surface). One can plot the mean/max/min of the z value in a 2D heatmap plot
with 20*20 bins with

    xyzPlots(x, y, z, measure='mean', bins = 20)

see example.py that generates

![](/plots/example.pdf)
