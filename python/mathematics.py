## @package mathematics
#  Library for rotation matrices and other pure math function used in aerospace engineering
#  
#  @author Griffin Jourda
import numpy as np

## rot_1
#  @brief rotation about x-axis
# 
#  This function computes a rotation matrix around the x-axis. The function uses the passive counterclockwise 
#  convention
#  @param theta rotation angle (rad) 
#  @return dcm rotation matrix
def rot_1(theta):
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[1, 0, 0],
                    [0, c, s],
                    [0, -s, c]])

    return dcm

## rot_2
#  @brief rotation about y-axis
# 
#  This function computes a rotation matrix around the y-axis. The function uses the passive counterclockwise 
#  convention
#  @param theta rotation angle (rad) 
#  @return dcm rotation matrix
def rot_2(theta):
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[c, 0, -s],
                   [0, 1, 0],
                   [s, 0, c]])

    return dcm

## rot_3
#  @brief rotation about z-axis
# 
#  This function computes a rotation matrix around the z-axis. The function uses the passive counterclockwise 
#  convention
#  @param theta rotation angle (rad) 
#  @return dcm rotation matrix
def rot_3(theta):
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[c, s, 0],
                   [-s, c, 0],
                   [0, 0, 1]])

    return dcm

## skew_sym
#  @brief create skew-symmetric matrix from vector 
# 
#  This function computes the skew-symmetric "cross-product" matrix of a vector, such that skew_sym(v1) @ v2 == cross(v1, v2)
#  @param vec: vector to compute skew-symmetric matrix of 
#  @return skew: skew-symmetric matrix of input vector
def skew_sym(vec):
    skew = np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])

    return skew