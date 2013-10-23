#! /usr/bin/python

import numpy as np
import sys
import math
import scipy.interpolate as interp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from processing import *
#Load in arguments
from arg_loader import *

#define functions used in rest of program

#dot product function
def dot_product(coord_a,coord_b):
	if len(coord_a) != len(coord_b):
		print "Cannot do dot product, program exiting..."
		sys.exit()
	dot_product = 0.0
	for i in range(len(coord_a)):
		dot_product+=coord_a[i]*coord_b[i]
	return dot_product

#Contour generation functions
def profile(distance_along_lim_contour): 
#TODO: find where the zero point is!
	return [(- 0.03112*x - 4.62E-04*x**2 + 1.46E-06*x**3 -1.15E-09*x**4) for x in distance_along_lim_contour]

def mass_change(mass): #uses location of hump
	return 11628.59959/((mass-550)**0.91498)
	
#CLS functions
def m12_cls(m12):
	return (2.5463622E-14*m12**5 - 1.5308754E-010*m12**4 + 3.6366585E-07*m12**3 - 4.2735866E-04*m12**2 + 2.4873059E-01*m12 - 5.6413813E+01)
	
def mg_cls(mg):
	return (1.826E-10*mg**3-6.511E-7*mg**2+6.2152E-4*mg+0.83034)

#Which contour to draw (in percentage)
contour_lim = sys.argv[9]
print "Limit contour calculated is "+contour_lim
#Which input file contains the group of masses we are interested in
input_file = sys.argv[10]
print "Input file used is "+input_file

#load in masses
input_masses = np.loadtxt(input_file,delimiter=',')

#Generate m3 master line
path_temp=processing(master_file,x_u,y_u,x_l,y_l,3000,x_new_res,y_new_res,coord_opt,int(contour_lim))
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
	
path_master_tck=interp.splrep(path_master[:,0],path_master[:,1],k=2,s=0)#no smoothing is needed as points are already extremely close
path_master_deriv=interp.splev(path_master[:,0],path_master_tck,der=1)

#distance from topmost point on master 
#TODO: Make this more general
distance_along_master=[0]
for i in range(1,len(path_master)):
	distance_along_master.append(distance_along_master[i-1]+np.sqrt((path_master[i][0]-path_master[i-1][0])**2+(path_master[i][1]-path_master[i-1][1])**2))
distance_profile = profile(distance_along_master)

#Dictionary for m3 masses and contours
m3_dict = {}

#Generate m3 lines
for m3_mass in input_masses[:,2]:
	#if contour for m3_mass already drawn, do not draw it again
	if m3_mass in m3_dict:
		continue
	overall_shift = mass_change(m3_mass)
	new_contour_dist = [ x+overall_shift for x in distance_profile ] 
	new_contour=[]
	for i in range(len(path_master)):
		y=new_contour_dist[i]*math.sin(math.atan(-1.0/path_master_deriv[i]))
		x=new_contour_dist[i]*math.cos(math.atan(-1.0/path_master_deriv[i]))
		new_contour.append([path_master[i][0]+x, path_master[i][1]+y])
	m3_dict[m3_mass]=new_contour

#Find point on m3 curve at which the normal to (m12,mg) is drawn
for masses in input_masses:
	#Load in values for m12 and mg, and look up the contour for m3
	m12 = masses[1]
	mg = masses[0]
	m3_contour = m3_dict[masses[2]]
	
	#Calculate tck values and first derivative for m3 contour
	m3_contour_tck=interp.splrep([x[0] for x in m3_contour],[x[1] for x in m3_contour],k=2,s=0)#no smoothing is needed as points are already extremely close
	m3_contour_deriv=interp.splev([x[0] for x in m3_contour],m3_contour_tck,der=1)
	
	#Calculate vector joining (m12,mg) to first point on m3 contour
	join_vec = [m12-m3_contour[0][0],mg-m3_contour[0][1]]
	#Calcuate vector of tangent at first point on m3 contour
	tangent_vec = [1,1*m3_contour_deriv[0]]
	#Calculate cosine of angle between vector joining (m12,mg) to first point on m3 contour to the vector of the tangent there
	cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
	#Initialise the sign of the cosine
	sign = cos/abs(cos)
	
	#Repeat above for all points on the m3 contour, until a sign change in the cosine is found
	#This indicates a point where the cosine = 0 i.e. where the vector between (m12,mg) and a point on the m3 contour and the tangent there are 
	#perpendicular. At this point the normal connecting the contour to (m12,mg) has been found.
	#The average of the points about the sign change is taken as an estimate of where the cosine=0 point is.
	for i in range(1,len(m3_contour)):
		join_vec = [m12-m3_contour[i][0],mg-m3_contour[i][1]]
		tangent_vec = [1,1*m3_contour_deriv[i]]
		cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
		sign_temp = cos/abs(cos)
		if sign_temp != sign:
			break
		sign = sign_temp
	if i == len(m3_contour)-1 and sign == sign_temp:
		print "No CLS calculated because the program was not able to draw a normal from the m3 curve to the point"
	else:
		x = (m3_contour[i][0]+m3_contour[i-1][0])/2
		y = (m3_contour[i][1]+m3_contour[i-1][1])/2
	print x
	print y
	
	cls_x = m12_cls(m12)/m12_cls(x)
	cls_y = mg_cls(mg)/mg_cls(y)
	CLS = cls_x*cls_y
	if CLS*100 > 100:
		CLS = 100.0
	print CLS

#Sanity check plots
	#gradient = (mg-y)/(m12-x)
	#X = np.mgrid[x_l:x_u:x_new_res*1j]
	#Y= [ element*gradient+(mg-gradient*m12) for element in X ] 
	#plt.ylim(y_l,y_u)
	#plt.xlim(x_l,x_u)
	#plt.plot(X,Y)
	#plt.plot(m12,mg,'ro')
	#plt.plot([x[0] for x in m3_contour], [x[1] for x in m3_contour])
	#plt.show()
