#! /usr/bin/python

import numpy as np
import sys
import math
import scipy.interpolate as interp
import matplotlib.pyplot as plt

#Use function processing
from processing import *
#Loads in arguments
from arg_loader import *

def profile(distance_along_lim_contour):
	return [(-43.55077 + 0.2325*x -8.88E-04*x**2 + 2.35E-06*x**3 + -2.34E-09*x**4) for x in distance_along_lim_contour]
	
def mass_change(mass): #uses location of hump
	return -320.05879 + 5.07453*mass -0.01201*mass**2 + 1.02E-05*mass**3 -2.93E-09*mass**4


mass = float(sys.argv[9]) #mass of particle for contour that needs to be drawn

#Analysis for limit contour (master line)
path_temp=processing(master_file,x_u,y_u,x_l,y_l,3000,x_new_res,y_new_res,coord_opt)
if isinstance(path_temp,np.ndarray)==False:
	if path_temp == -1:
		print "Master file not found. Program exiting..."
		sys.exit()
	if path_temp == -2:
		print "Master file has no data that reaches 95%. Program exiting..."
		sys.exit()
	if path_temp == -3:
		print "Master file data co-ordinates not found. Program exiting..."
		sys.exit()
else:
	path_master = path_temp
	
#TODO: still undertake interpolation despite duplicated points in x co-ordinates
#TODO: GET FUNCTIONAL FORM FOR MASTER PATH! (FUNCTIONAL SURFACE FOR SPLINE FROM COEFFS, THEN CROSS-SECTION /W SURFACE)
path_master_tck=interp.splrep(path_master[:,0],path_master[:,1],k=2,s=0)#no smoothing is needed as points are already extremely close
path_master_deriv=interp.splev(path_master[:,0],path_master_tck,der=1)

#distance from topmost point on master. As we do not have a functional form, distances are added linear piecewise.
distance_along_master=[0]
for i in range(1,len(path_master)):
	distance_along_master.append(distance_along_master[i-1]+np.sqrt((path_master[i][0]-path_master[i-1][0])**2+(path_master[i][1]-path_master[i-1][1])**2))

distance_profile = profile(distance_along_master)
overall_shift = mass_change(mass)
new_contour_dist = [ x+overall_shift for x in distance_profile ] 

new_contour=[]

for i in range(len(path_master)):
		y=new_contour_dist[i]*math.sin(math.atan(-1.0/path_master_deriv[i]))
		x=new_contour_dist[i]*math.cos(math.atan(-1.0/path_master_deriv[i]))
		new_contour.append([path_master[i][0]+x, path_master[i][1]+y])

plt.plot(np.array(new_contour)[:,0],np.array(new_contour)[:,1])
plt.show()






















