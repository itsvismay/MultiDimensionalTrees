from __future__ import division
import os
import sys
import re
from math import *

#Define ray tracing resolution
ray_x=640 #Use 1280 for high quality or 320 for testing
ray_y=512 #Use 1024 for high quality or 256 for testing

#Define ray trace options
ray_mode=0 
ray_shadows=0
ray_antialias=1 #Set to 0 for faster processing
ray_ambient=0.3 #Comment out with '#' sign here and in ray() to use the default light 
ray_save_dpi=72

#Start frame counter
frame_number=1

#Define the method to ray trace the frame and dump image to the hard drive
def ray():
    global ray_x, ray_y, ray_mode, ray_antialias, ray_save_dpi, frame_number
    #Set ray tracer options
    cmd.set("ray_shadows", ray_shadows)
    cmd.set("ambient", ray_ambient)#Comment out with '#' sign to use the default light
    cmd.set("ray_trace_mode", ray_mode)
    #Run ray tracer
    cmd.ray(ray_x, ray_y, ray_antialias)
    #Save image
    cmd.png("mov"+`frame_number`, dpi=ray_save_dpi)
    frame_number+=1

#Set the camera view accoring to 18-point matrix
#Use get_view method in PyMol command line to get the matrix
#@x 18-point PyMol view matrix
def setView (x):    
    view=""
    for i in range(len(x)):
        view+=`x[i]`
        if (i != (len(x)-1)):
            view+=","
    cmd.set_view(view)
    
#Rotate an object 360 degrees
#@frames: number of frames to animate (integer)
#@rotation_axis: axis to rotate the structure ('x', 'y' or 'z')
#@originpoint: center of rotation 
#(can be a PyMol object or a particular residue/atom)
def roll (frames, rotation_axis, originpoint):
    #Set rotation center to user-defined "originpoint"
    #cmd.origin(originpoint)
    #Define the angle to rotate per frame
    angle=360/frames
    frame=1
    #Loop over frames, make a small rotation, 
    #ray trace the scene and save the image
    while frame <= frames:
        #Use PyMol rotation function
        cmd.rotate(rotation_axis, angle)
        #Ray trace and save current scene
        ray()
        frame+=1
roll(100, 'x',(30,0,0))
