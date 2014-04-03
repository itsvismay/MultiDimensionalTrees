#PyMol script to animate structural alignment of two biological macromolecules
#It contains a set of basic methods and an example of a user-defined scenario 
#Copyright (C) 2013 Dmitry Suplatov d.a.suplatov@belozersky.msu.ru http://biokinet.belozersky.msu.ru

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import os
import sys
import re
from math import *

#Define ray tracing resolution
ray_x=640 #Use 2560 for high quality or 320 for testing
ray_y=512 #Use 2048 for high quality or 256 for testing

#Define ray trace options
ray_mode=0 #Use 0 for "regular" type images
ray_shadows=0
ray_antialias=1 #Set to 0 for faster processing
ray_ambient=0.3 #Comment out with '#' sign here and in ray() to use the default light
ray_save_dpi=72

#Start frame counter
frame_number=1

def ray():
    global ray_x, ray_y, ray_mode, ray_antialias, ray_save_dpi, frame_number
    cmd.set("ray_shadows", ray_shadows)
    cmd.set("ambient", ray_ambient) #Comment out with '#' sign to use the default light 
    cmd.set("ray_trace_mode", ray_mode)
    cmd.ray(ray_x, ray_y, ray_antialias)
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
    cmd.origin(originpoint)
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

#Rotate two different objects 360 degrees aroung thier axes
#@frames: number of frames to animate (integer)
#@rotation_axis: axis to rotate the structure ('x', 'y' or 'z')
#@obj1: first object (protein molecule)
#@obj2: second object (protein molecule)
def two_roll (frames, rotation_axis, obj1, obj2):
    angle=360/frames        
    frame=1
    #Loop over frames, make small rotation to each object, 
    #ray trace the scene and save the image
    while frame <= frames:
        #Set rotation center to object1
        cmd.origin(obj1)
        #Rotate object1
        cmd.rotate(rotation_axis, angle, obj1)
        #Set rotation center to object2
        cmd.origin(obj2)
        #Rotate object2
        cmd.rotate(rotation_axis, angle, obj2)
        #Ray trace and save the image
        ray()
        frame+=1

#Simulates alignment of two protein structures.
#The method should be used on two already superimposed molecules
#@frames_roll - number of frames to animate rotation of two de-alinged molecules
#@frames_align - number of frames to animate superimposition
#@obj1: first object (protein molecule)
#@shift_x1: shift the first object to the left to de-align
#@obj2: second object (protein molecule)
#@shift_x2: shift the second object to the right to de-align
#@rotation_axis: axis to rotate the structures during animation ('x', 'y' or 'z')
def align(frames_roll, frames_align, obj1, shift_x1, obj2, shift_x2, rotation_axis):    
    #De-align protein molecules
    #Shift to the left with "-" sign
    cmd.translate("\"[-"+`shift_x1`+", 0, 0]\"", obj1)
    #Shift to the right
    cmd.translate("\"["+`shift_x2`+", 0, 0]\"", obj2)
    #Take the first snapshot of the molecules in de-aligned state
    ray()    

    #Roll two de-aligned molecules in place
    two_roll(frames_roll, 'y', '1tcb', '1j2e')

    #Animate superimposition
    angle=360/frames_align
    frame=1
    while (frame <= frames_align):
        #Define translation per frame
        cmd.translate("\"["+`shift_x1/frames_align`+", 0, 0]\"", obj1)
        cmd.translate("\"[-"+`shift_x2/frames_align`+", 0, 0]\"", obj2)
        
        #Define rotation per frame
        cmd.origin(obj1)
        cmd.rotate(rotation_axis, angle, obj1)
        cmd.origin(obj2)
        cmd.rotate(rotation_axis, angle, obj2)
        
        #Create the image
        ray()
        frame+=1

#Changes position of the camera
#@frames: number of frames to animate
#@x: view matrix of the starting frame
#@y: view matrix of the target frame
def changeView (frames,x,y):        
    
    setView(x)
    ray()
    
    #get quaternions for interpolation between starting and target scenes
    qstart = mat2quat(x[0:9])
    qend = mat2quat(y[0:9])

    frame=1
    while frame <= frames:                    
        n=[]
        
        qcur = slerp(qstart, qend, frame/frames)
        matcur = quat2mat(qcur)

        for i in range(9):
            update=matcur[i]
            if abs(update) < 0.001:
                update = 0
            n.insert(i, update)
	
	for i in range(len(x)):	    
            if (i>8):
    	        update=x[i] + (y[i] - x[i]) * (frame/frames)
    		if abs(update) < 0.001:
            	    update = 0
                n.insert(i, update)

        setView(n)
        frame+=1    
        ray()
 
#-------------------------------------------------------------------------------       
#Small library of math methods to perform Spherical Linear Interpolation (SLERP)
#The code has been translated to python from sources available at
#http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
#-------------------------------------------------------------------------------
def slerp(qa, qb, t):
    qm=[]
    
    #Calculate angles between quaternions
    cosHalfTheta = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3]
    #if qa=qb or qa=-qb then theta = 0 and we can return qa
    
    if (cosHalfTheta < 0):
        for i in range(4):
            qb[i] = -qb[i];
        cosHalfTheta = -cosHalfTheta
    
    if (abs(cosHalfTheta) >= 1.0):        
        for i in range(4):
            qm.insert(i,qa[i])
        return qm
    
    #Calculate temporary values
    halfTheta = acos(cosHalfTheta)
    sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)
    
    if (fabs(sinHalfTheta) < 0.000005):
    #fabs is floating point absolute
        for i in range(4):
            qm.insert(i, qa[i] * 0.5 + qb[i] * 0.5)        
        return qm

    ratioA = sin((1 - t) * halfTheta) / sinHalfTheta
    ratioB = sin(t * halfTheta) / sinHalfTheta
    #calculate Quaternion.
    for i in range(4):
        qm.insert(i, qa[i] * ratioA + qb[i] * ratioB)
    return qm

def mat2quat(m):

    tr = m[0] + m[4] + m[8]
    
#    m00=0 m01=1 m02=2
#    m10=3 m11=4 m12=5
#    m20=6 m21=7 m22=8
    
    if (tr > 0):
      S = sqrt(tr+1) * 2;
      qw = 0.25 * S;
      qx = (m[7] - m[5]) / S;
      qy = (m[2] - m[6]) / S; 
      qz = (m[3] - m[1]) / S; 
    
    elif ((m[0] > m[4])&(m[0] > m[8])):
      S = sqrt(1 + m[0] - m[4] - m[8]) * 2
      qw = (m[7] - m[5]) / S;
      qx = 0.25 * S;
      qy = (m[1] + m[3]) / S; 
      qz = (m[2] + m[6]) / S; 
    
    elif (m[4] > m[8]):
      S = sqrt(1 + m[4] - m[0] - m[8]) * 2
      qw = (m[2] - m[6]) / S;
      qx = (m[1] + m[3]) / S; 
      qy = 0.25 * S;
      qz = (m[5] + m[7]) / S; 
    
    else:
      S = sqrt(1 + m[8] - m[0] - m[4]) * 2
      qw = (m[3] - m[1]) / S;
      qx = (m[2] + m[6]) / S;
      qy = (m[5] + m[7]) / S;
      qz = 0.25 * S;
    
    return [qx,qy,qz,qw]
        

def quat2mat( Q ):
    
    xx = Q[0]*Q[0]
    xy = Q[0]*Q[1]
    xz = Q[0]*Q[2]
    xw = Q[0]*Q[3]
    yy = Q[1]*Q[1]
    yz = Q[1]*Q[2]
    yw = Q[1]*Q[3]
    zz = Q[2]*Q[2]
    zw = Q[2]*Q[3]
    
    
    M = [1.0 - 2*yy - 2*zz,
         2*xy - 2*zw,
         2*xz + 2*yw,
        2*xy + 2*zw,
        1 - 2*xx - 2*zz,
        2*yz - 2*xw,
        2*xz - 2*yw,
        2*yz + 2*xw,
        1 - 2*xx - 2*yy]
    
    return M

#---------------------------------
#USER-DEFINED SCENARIO STARTS HERE
#---------------------------------

#Set starting viewpoint
#Type get_view in PyMol command line and copy-paste the 18-point matrix.
#The view should be tested in PyMol viewer and set accordingly to capture the important part of the scene
setView([0.913405001,    0.009108904,   -0.404045910,
     0.199869096,    0.857655287,    0.471460521,
     0.351102561,   -0.511840045,    0.782457471,
    -0.186222672,   -1.656485915, -355.204223633,
    45.377914429,   22.392894745,   81.392692566,
   300.970520020,  402.013854980,  -20.000000000])

#Animate structural alignment
#The method starts from two aligned objects
#Then, the two objects - 1tcb and 1j2e - are each moved 30 points in different directions. First frame is captured.
#   Be sure to set the view (setView command above) so that it covers the scene
#The two objects start rolling in place for 90 frames
#The two objects continue rolling but start to move to each other. This "alignment" step takes another 45 frames
align(90, 45, '1tcb', 30, '1j2e', 30, 'y')

#Since the two objects are now aligned, we need to zoom in as we dont want to much white on the edges of our movie
#The transition will be make in 30 frames (first number passed to the method)
changeView(30,
    [0.913405001,    0.009108904,   -0.404045910,
     0.199869096,    0.857655287,    0.471460521,
     0.351102561,   -0.511840045,    0.782457471,
    -0.186222672,   -1.656485915, -355.204223633,
    45.377914429,   22.392894745,   81.392692566,
   300.970520020,  402.013854980,  -20.000000000],
    
    [0.913405001,    0.009108904,   -0.404045910,
    0.199869096,    0.857655287,    0.471460521,
     0.351102561,   -0.511840045,    0.782457471,
     0.000000000,    0.000000000, -256.684509277,
    28.279685974,   22.153846741,   77.451515198,
   206.162841797,  307.206085205,  -20.000000000])    

#Now the two objects are aligned and zoomed in - time to make, yet again, a 360-degree rotation
roll(120, 'y', '1tcb')

#Zoom on the active site
#The first 18-point matrix represents the coordinated of the last scene
#The second 18-point matrix was selected manually by observing the alignment in PyMol viewer,
#   selecting appropriate scene and using get_view command
changeView(45,
     [0.913405001,    0.009108904,   -0.404045910,
     0.199869096,    0.857655287,    0.471460521,
     0.351102561,   -0.511840045,    0.782457471,
     0.000000000,    0.000000000, -256.684509277,
    28.279685974,   22.153846741,   77.451515198,
   206.162841797,  307.206085205,  -20.000000000],
     
    [-0.245052531,    0.443468809,    0.860734820,
     0.717557430,    0.679393291,   -0.145888939,
    -0.650069058,    0.582493782,   -0.485329747,
     0.000000000,    0.000000000,  -44.810520172,
    34.562198639,   28.321912766,   79.162284851,
    39.730285645,   49.890754700,  -19.999998093])

#Final rotation to explore the active site.
#NE2 atom of the catalytic His224 will be used as rotation center around the Y-axis
roll(120, 'y', '1tcb and resid 224 and name NE2')
