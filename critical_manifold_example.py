# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 11:50:39 2020

Minimal code required to produce critical manifold curves.

"""

# Call necessary libraries. lubrication.py is available on our repository https://github.com/youngmp/park_fai_2020. Numpy and matplotlib are available on https://www.scipy.org/
import numpy as np
import matplotlib.pyplot as plt
from lubrication import lubrication

# Generate a figure object and create a subplot which we will populate with the critical manifold
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)

# Generate an instance of the lubrication model using parameters from Fig. 2.2A
a = lubrication(phi1=.57,Rp=.96,Rc=1.22,
                pi3=1,pi4=4.7,pi5=0.1,pi6=10,
                mu=1.2,T=3,constriction='piecewise',U0=0.2,
                dt=0.02,eps=1,
                F0=50)

# Define the number of grid points in each axis of the U-Z domain. The total number of grid points on the 2D mesh will be M*M.
M = 200

# Define the velocity and position ranges
U = np.linspace(-.4,.4,M)/a.U_scale # velocity range
Z = np.linspace(-5,5,M)/a.Z_scale # position range

# Create the meshgrid object for 3D or contour plotting. See https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html
Zm,Um = np.meshgrid(Z,U)

# The viscous drag function in the lubrication model object only accepts scalar or 1D-array inputs, so we extract the values of zeta for each position gridpoint
ze = np.zeros((M,M))
for i in range(M):
    for j in range(M):
        ze[i,j] = a.viscous_drag(Zm[i,j])

# Define the surface from which we will extract the zero contour
C_0 = a.F(Um)-ze*Um

# Draw the zero contour lines of the surface C_0 using the 'levels' option '[0]'. See https://matplotlib.org/api/_as_gen/matplotlib.pyplot.contour.html for more details.
# a.Z_scale and a.U_scale redimensionalize the axes.
# ax_contour saves information about the contour plot. We will extract the contour U-Z positions from this variable.
ax_contour = ax.contour(Zm*a.Z_scale,Um*a.U_scale,C_0,[0])

# Clear the uncolored, unordered contour lines
ax.clear()

######################
# To color the contour lines as stable and unstable, much manual work must be done. First, we note beforehand, based on the bifurcation diagrams, which manifolds are stable are unstable. We then use Python index slicing to extract the corresponding data points from ax_contour. Separating the data points is necessary in order to plot the points as different colors. While color with a scatter plot would be easier, they are poorly optimized in Python and tend to make PDFs load slowly.
######################

# Extract all contour lines (there only happen to be 3 in our figure).
p0=ax_contour.collections[0].get_paths()[0]
p1=ax_contour.collections[0].get_paths()[1]
p2=ax_contour.collections[0].get_paths()[2]

# By plotting the different contour lines in a separate program, we find that i=0 (p0) is the right-most curve, i=1 (p1) is the left-most curve, and i=2 (p2) is the upper curve.

# Extract the vertices of the contour lines for plotting.
x0=p0.vertices[:,0];y0=p0.vertices[:,1]
x1=p1.vertices[:,0];y1=p1.vertices[:,1]
x2=p2.vertices[:,0];y2=p2.vertices[:,1]

# From the bifurcation figures, we note that contours i=0 and i=1 have folds where they are stable below the fold and unstable above. So we split the curve into above and below parts so we can plot and color them differently.

# Find the x-coordinate max of p1 and the x-coordinate min of p0 and use the index to mark the y-coordinate value of the fold.
idx0 = np.argmin(x0)
idx1 = np.argmax(x1)

# Extract the stable and unstable coordinates of the left and right contours.
x_right_s,y_right_s = x0[:idx0],y0[:idx0] # right stable branch
x_right_u,y_right_u = x0[idx0:],y0[idx0:] # right unstable branch    

x_left_s,y_left_s = x1[idx1:],y1[idx1:] # left stable branch
x_left_u,y_left_u = x1[:idx1],y1[:idx1] # left unstable branch

# Plot the contours with colors corresponding to stability
ax.plot(x_right_s,y_right_s,color='k',label='Stable Manifold',lw=2)
ax.plot(x_right_u,y_right_u,color='tab:red',label='Unstable Manifold',lw=2)

ax.plot(x_left_s,y_left_s,color='k',lw=2)
ax.plot(x_left_u,y_left_u,color='tab:red',lw=2)

ax.plot(x2,y2,color='k',lw=2)

ax.legend()

plt.savefig('critical_manifold_example.png')