#! /usr/bin/python

import numpy as np
import sys
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import re
from processing import *
#Load in arguments
from arg_loader import *

def m3_cls(m3):
	return (-5.08296E-12*m3**4+1.8E-8*m3**3-2.3065E-5*m3**2+0.01253*m3-1.4452)

def m12_cls(m12):
	return (2.5463622E-14*m12**5 - 1.5308754E-010*m12**4 + 3.6366585E-07*m12**3 - 4.2735866E-04*m12**2 + 2.4873059E-01*m12 - 5.6413813E+01)
	
def mg_cls(mg):
	return (1.826E-10*mg**3-6.511E-7*mg**2+6.2152E-4*mg+0.83034)

def dot_product(coord_a,coord_b):
	if len(coord_a) != len(coord_b):
		print "Cannot do dot product, program exiting..."
		sys.exit()
	dot_product = 0.0
	for i in range(len(coord_a)):
		dot_product+=coord_a[i]*coord_b[i]
	return dot_product

print master_file
values_orig = np.genfromtxt(master_file,delimiter=',')
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
path_master_tck=interp.splrep(path_master[:,0],path_master[:,1],k=2,s=0)#no smoothing is needed as points are already extremely close
path_master_deriv=interp.splev(path_master[:,0],path_master_tck,der=1)

#Set up output file
#Get mass of particle from master data file
mass = (re.search('(?<=/)\d+',master_file)).group(0)
#Get type of particle from master data file
try:
	master_type=(re.search('m\w+',master_file)).group(0)
except AttributeError:
	master_type=(re.search('m\d+',master_file)).group(0)

outfile = 'plane_scan_'+master_type+'_'+mass+'.txt'	
fout = open(outfile,'w')
if master_type == 'mg':
	print "Rows are m3, columns are m12"
	x_label = 2
	y_label = 1
elif master_type == 'm12':
	print "Rows are m3, columns are mg"
	x_label = 2
	y_label = 0
elif master_type == 'm3':
	print "Rows are m12, columns are mg"
	x_label = 1
	y_label = 0
else:
	print "No acceptable master data type found. Program exiting..."
	sys.exit()
	
#Load in co-ordinates for x and y values
if coord_opt == 'n':
	try:
   		with open(master_file+"_coord"): pass
	except IOError:
  			print "Co-ordinates for master file not found. Program exiting..."
  			sys.exit()
	else:
		data_coord=open(master_file+"_coord",'r')
		x_mass_matrix=((data_coord.readline()).strip()).split(',')
		x_mass_matrix = [float(i) for i in x_mass_matrix ]
		y_mass_matrix=((data_coord.readline()).strip()).split(',')
		y_mass_matrix= [float(i) for i in y_mass_matrix ]
elif coord_opt == 'd':
	x_mass_matrix = np.mgrid[x_l:x_u:len(values_orig[0])*1j]
	y_mass_matrix = np.mgrid[y_l:y_u:len(values_orig)*1j]

#Make list of functions needed for cls
cls_func=[mg_cls,m12_cls,m3_cls]

#Find tangent between co-ordinates and points on path_master
for y_mass in y_mass_matrix:
	for x_mass in x_mass_matrix:
		#Calculate angle between tangent and line joining two points on the path
		join_vec = [x_mass-path_master[0][0],y_mass-path_master[0][1]]
		tangent_vec = [1,1*path_master_deriv[0]]
		cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
		sign = cos/abs(cos)
		for i in range(len(path_master)):
			join_vec = [x_mass-path_master[i][0],y_mass-path_master[i][1]]
			tangent_vec = [1,1*path_master_deriv[i]]
			cos = dot_product(join_vec,tangent_vec)/(np.sqrt(dot_product(join_vec,join_vec))*np.sqrt(dot_product(tangent_vec,tangent_vec)))
			sign_temp = cos/abs(cos)
			if sign_temp != sign:
				break
			sign = sign_temp

		if i == len(path_master)-1 and sign == sign_temp:
			if x_mass_matrix.index(x_mass) == len(x_mass_matrix)-1:
				fout.write('nan')
			else:
				fout.write('nan,')
			continue
		else:
			x = (path_master[i][0]+path_master[i-1][0])/2
			y = (path_master[i][1]+path_master[i-1][1])/2
			#Plotting the normal...pick points; this is essentially a test
			if x_mass == 950 and y_mass == 1000:
				gradient = (y_mass-y)/(x_mass-x)
				X = np.mgrid[x_l:x_u:x_new_res*1j]
				Y= [ element*gradient+(y_mass-gradient*x_mass) for element in X ] 
				plt.ylim(y_l,y_u)
				#plt.plot(X,Y)
				#plt.plot(x_mass,y_mass,'ro')
		#Choose expression
			cls_x = cls_func[x_label](x_mass)/cls_func[x_label](x)
			cls_y = cls_func[y_label](y_mass)/cls_func[y_label](y)
			CLS = cls_x*cls_y
			if CLS*95 > 100:
				CLS = 100.0/95.0
			if x_mass_matrix.index(x_mass) == len(x_mass_matrix)-1:
				fout.write(str(CLS*95))
			else:
				fout.write(str(CLS*95)+',')
	fout.write('\n')
	
fout.close()
values_gen = np.genfromtxt(outfile,delimiter=',')
#Change this line depending on what you want to see...
plt.clf()
v=np.linspace(90,99,3,endpoint=True)
fig = plt.figure()
gen_plot = fig.add_subplot(1,2,1)
plt.contour(x_mass_matrix,y_mass_matrix,values_gen,v)
orig_plot = fig.add_subplot(1,2,2)
plt.contour(x_mass_matrix,y_mass_matrix,values_orig,v)
plt.colorbar()
plt.show()
