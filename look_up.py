#! /usr/bin/python

import numpy as np
import sys
import math
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from processing import *
#Load in arguments
from arg_loader import *

#dot product function
def dot_product(coord_a,coord_b):
	if len(coord_a) != len(coord_b):
		print "Cannot do dot product, program exiting..."
		sys.exit()
	dot_product = 0.0
	for i in range(len(coord_a)):
		dot_product+=coord_a[i]*coord_b[i]
	return dot_product

#load in masses

mg = float(sys.argv[9])
m12 = float(sys.argv[10])
m3  = float(sys. argv[11])

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

#Generate m3 master line
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

#Generate m3 line
path_master_tck=interp.splrep(path_master[:,0],path_master[:,1],k=2,s=0)#no smoothing is needed as points are already extremely close
path_master_deriv=interp.splev(path_master[:,0],path_master_tck,der=1)

#distance from topmost point on master #TODO: Make this more general
distance_along_master=[0]
for i in range(1,len(path_master)):
	distance_along_master.append(distance_along_master[i-1]+np.sqrt((path_master[i][0]-path_master[i-1][0])**2+(path_master[i][1]-path_master[i-1][1])**2))

distance_profile = profile(distance_along_master)
overall_shift = mass_change(m3)
new_contour_dist = [ x+overall_shift for x in distance_profile ] 

new_contour=[]

for i in range(len(path_master)):
		y=new_contour_dist[i]*math.sin(math.atan(-1.0/path_master_deriv[i]))
		x=new_contour_dist[i]*math.cos(math.atan(-1.0/path_master_deriv[i]))
		new_contour.append([path_master[i][0]+x, path_master[i][1]+y])

new_contour_tck=interp.splrep([x[0] for x in new_contour],[x[1] for x in new_contour],k=2,s=0)#no smoothing is needed as points are already extremely close
new_contour_deriv=interp.splev([x[0] for x in new_contour],new_contour_tck,der=1)
		
#Calculate CLS

#Find point on m3 curve at which the normal to (m12,mg) is drawn
join_vec = [m12-new_contour[0][0],m3-new_contour[0][1]]
tangent_vec = [1,1*new_contour_deriv[0]]
cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
sign = cos/abs(cos)
for i in range(len(path_master)):
	join_vec = [m12-new_contour[i][0],mg-new_contour[i][1]]
	tangent_vec = [1,1*new_contour_deriv[i]]
	cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
	sign_temp = cos/abs(cos)
	if sign_temp != sign:
		break
	sign = sign_temp
if i == len(new_contour)-1 and sign == sign_temp:
	print "No CLS calculated because the program was not able to draw a normal from the m3 curve to the point"
else:
	x = (new_contour[i][0]+new_contour[i-1][0])/2
	y = (new_contour[i][1]+new_contour[i-1][1])/2

print x
print y

cls_x = m12_cls(m12)/m12_cls(x)
cls_y = mg_cls(mg)/mg_cls(y)
CLS = cls_x*cls_y
if CLS*95 > 100:
	CLS = 100.0/95.0

print CLS
plt.plot([x[0] for x in new_contour],[x[1] for x in new_contour])
plt.show()
