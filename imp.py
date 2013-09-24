#! /usr/bin/python

import numpy as np
import sys
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import re

#Use function processing
from processing import *
#Loads in arguments
from arg_loader import *

print master_file

output = open("output.txt",'w')

def sign_finder(num):
	return num/abs(num)
	
#Analysis for master line
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

#Set up for 3D plot
#y
mass_coord=[]
#x
#distance from topmost point on master. As we do not have a functional form, distances are added linear piecewise.
distance_along_master=[0]
for i in range(1,len(path_master)):
	distance_along_master.append(distance_along_master[i-1]+np.sqrt((path_master[i][0]-path_master[i-1][0])**2+(path_master[i][1]-path_master[i-1][1])**2))	
#z
distance_coord=[]

#TODO: plot contour /w even spread x data points

#Analysis of data files provided
for filename in sys.argv[9:]:
	data_temp=processing(filename,x_u,y_u,x_l,y_l,3000,x_new_res,y_new_res,coord_opt)
	if isinstance(data_temp,np.ndarray)==False:
		if data_temp == -1:
			print "Data file "+filename+" not found. Skipping..."
			continue
		if data_temp == -2:
			print "Data file "+filename+" has no data reaching 95%. Skipping..."
			continue
		if data_temp == -3:
			print "Data file "+filename+" co-ordinates not found. Skipping..."
			continue
	else:
		path_data=data_temp
		
	#find mass from data file
	#mass = filename.split('_')[0]
	mass = (re.search('(?<=/)\d+', filename)).group(0)
	distance_coord_temp= []
	
	#get intersection of line with path_data
	for i in range(len(path_master)):
		norm_grad = -1.0/path_master_deriv[i]
		constant = path_master[i][1] - norm_grad*path_master[i][0]
		
	#Stores (possibly multiple) intersection points
		intersects_temp=[]
	#Initialise sign finding
		sign=sign_finder((path_data[0][1]-constant)/norm_grad-path_data[0][0])
		for j in range(len(path_data)):
			sign_temp=sign_finder((path_data[j][1]-constant)/norm_grad-path_data[j][0])
		#Look for changes of sign
			if (sign_temp != sign):
				intersects_temp.append(j)
		#update sign
			sign = sign_temp
		#If no intersection is found, write value of point as nan
		if len(intersects_temp) == 0:
			distance_intersect = np.nan
		else:
			distance_intersect = (x_u-x_l)*(y_u-y_l)
			for index in intersects_temp:
				x_intersect = (path_data[index-1][0]+path_data[index][0])/2
				y_intersect = (path_data[index-1][1]+path_data[index][1])/2
				distance_intersect_temp= np.sqrt((x_intersect-path_master[i][0])**2 + (y_intersect-path_master[i][1])**2)
			#Always use the shortest distance as the true intersection
				if distance_intersect_temp < distance_intersect:
					distance_intersect = distance_intersect_temp
					
	# Load in data for plot
		distance_coord_temp.append(distance_intersect)
	mass_coord.append(float(mass))
	distance_coord.append(distance_coord_temp)
	for i in range(len(distance_along_master)):
		output.write(str(distance_along_master[i])+'\t'+str(distance_coord_temp[i])+'\n')

#Make output plot
fig = plt.figure()

#PARAMETERIZE AS DISTANCE ALONG LIMIT CONTOUR
for i in range(len(mass_coord)):
	plt.plot(distance_along_master,distance_coord[i])
plt.show()
	




















