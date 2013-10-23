#! /usr/bin/python

#This file defines the processing function

import numpy as np
import scipy.interpolate as interp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def processing(filename,x_u,y_u,x_l,y_l,s_val,x_new_res,y_new_res,coord_opt,contour_lim):
	#load in data as a 2D matrix
	try:
   		with open(filename): pass
	except IOError:
  		return -1
	values = np.loadtxt(filename,delimiter=',')

	#Check if 95% limit will exist
	flag = False
	for row in values:
		for element in row:
			if element >= contour_lim:
				flag = True
				break
	if (flag == False):
		return -2
	
	#define data co-ordinates
	#TODO: take into account irregularly spaced data values
	if coord_opt == 'd':
		x = np.mgrid[x_l:x_u:len(values[0])*1j]
		y = np.mgrid[y_l:y_u:len(values)*1j]
	elif coord_opt == 'n':
	#request to read in co-ordinates noted in data file
		try:
   			with open(filename+"_coord"): pass
		except IOError:
  			return -3
  		else:
  			filename_coord = filename+"_coord"
  			data_coord=open(filename_coord)
  			x=((data_coord.readline()).strip()).split(',')
  			x = [float(i) for i in x ]
  			y=((data_coord.readline()).strip()).split(',')
  			y = [float(i) for i in y ]
	
	x,y = np.meshgrid(x,y)
	#interpolate using quadratic splines
	#Quadratic are used to better preserve asymptotic nature of plots
	#TODO:What value of s is optimal?
	tck = interp.bisplrep(x,y,values,kx=2,ky=2,s=s_val)
	
	#define points to interpolate over
	xnew,ynew = np.mgrid[x_l:x_u:(x_new_res*1j),y_l:y_u:(y_new_res*1j)]
	values_new = interp.bisplev(xnew[:,0],ynew[0,:],tck)

	#plot only the cls_level line
	v=np.linspace(contour_lim,contour_lim,2)
	cs = plt.contour(xnew,ynew,values_new,v)
	
	#Extract data of cls_level line
	#TODO: investigate syntax of this line
	#TODO: catch error where there is data below 95% but not enough to generate a contour
	return (cs.collections[0].get_paths()[0]).vertices
