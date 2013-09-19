#! /usr/bin/python

#This file contains a function that load in command line arguments
import sys

#1 master file
master_file = sys.argv[1]
print "Master file is "+master_file
#2 upper bound on x
x_u = float(sys.argv[2])
print "Upper bound of x is "+str(x_u)
#3 upper bound on y
y_u = float(sys.argv[3])
print "Upper bound of y is "+str(y_u)
#4 lower bound on x
x_l = float(sys.argv[4])
print "Lower bound of x is "+str(x_l)
#5 lower bound on y
y_l = float(sys.argv[5])
print "Lower bound of y is "+str(y_l)
#6 new resolution in x
x_new_res = float(sys.argv[6])
print "New resolution in x is "+str(x_new_res)
#7 new resolution in y
y_new_res = float(sys.argv[7])
print "New resolution in y is "+str(y_new_res)
#8 whether to generate co-ordinates or to read them from file
coord_opt = sys.argv[8]
print "Way to generate data co-ordinates is "+coord_opt
#9 onwards: data files to interpolate over
