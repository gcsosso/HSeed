#!/usr/bin/python

"""
(+) rotation matrix from one vector to another one
(+) rotation matrix to pole
(+) compute rotation matrix around given axis and angle
"""

import math
import numpy as np       # invert matrix

def rotation_matrix_v2v(vector, target):
    """ Roation matrix from one vector to another target vector.
 
    The solution is not unique as any additional rotation perpendicular to
    the target vector will also yield a solution)
     
    However, the output is deterministic.

    The algorithm works as following:
    (*) compute rotation matrix of vector to (1,0,0,...,n) (n - dimensional) ==> R1
    (*) compute rotation matrix of target to (1,0,0,...,n) (n - dimensional) ==> R2
    (*) the rotation R1.Transpose*R2 gives the roation matrix from vector to target

    """
 
    R1 = rotation_to_pole(target)
    R2 = rotation_to_pole(vector)
 
    return np.dot(R1.T, R2)
     
def rotation_to_pole(target):
    """ Rotate to 1,0,0,...,n
    
    works in n dimensions
    """
    # get dimension of vector
    n = len(target)
    working = target

    # make diagonal matrix of size n
    rm = np.eye(n)

    for i in range(1,n):
        angle = np.arctan2(working[0],working[i])
        rm = np.dot(rotation_matrix_inds(angle, n, 0, i), rm)
        working = np.dot(rm, target)
 
    return rm


def rotation_matrix(axis,theta):
    """
    compute rotation-matrix around 'axis' with angle 'theta'

    method: Euler-Rodriguez formula

    axis  --> numpy array
    theta --> in radians
    """

    theta = theta*(-1.0)
    axis = axis/math.sqrt(np.dot(axis,axis))
    a = math.cos(theta/2)
    b,c,d = -axis*math.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def rotation_matrix_inds(angle, n, ax1, ax2):
    """ 'n'-dimensional rotation matrix 'angle' radians in coordinate plane with
        indices 'ax1' and 'ax2' """
     
    s = np.sin(angle)
    c = np.cos(angle)
 
    i = np.eye(n)
 
    i[ax1,ax1] = s
    i[ax1,ax2] = c
    i[ax2,ax1] = c
    i[ax2,ax2] = -s
 
    return i
