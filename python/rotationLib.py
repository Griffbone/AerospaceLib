"""
Library for basic rotations and converting between attitude/orientation representations.

Author(s):
    Griffin Jourda

References:
    [1] https://www.vectornav.com/resources/inertial-navigation-primer/math-fundamentals/math-attitudetran
"""
import numpy as np


def rot_1(theta):
    """
    Return transformation matrix for a passive counterclockwise rotation about the x-axis.

    :param theta: angle of rotation (rad)
    :return dcm: rotation transformation matrix
    """
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[1, 0, 0],
                    [0, c, s],
                    [0, -s, c]])

    return dcm


def rot_2(theta):
    """
    Return transformation matrix for a passive counterclockwise rotation about the y-axis.

    :param theta: angle of rotation (rad)
    :return dcm: rotation transformation matrix
    """
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[c, 0, -s],
                   [0, 1, 0],
                   [s, 0, c]])

    return dcm


def rot_3(theta):
    """
    Return transformation matrix for a passive counterclockwise rotation about the z-axis.

    :param theta: angle of rotation (rad)
    :return dcm: rotation transformation matrix
    """
    c = np.cos(theta)
    s = np.sin(theta)

    dcm = np.array([[c, s, 0],
                   [-s, c, 0],
                   [0, 0, 1]])

    return dcm


def skew_sym(vec):
    """
    Return the skew-symmetric matrix of a vector.

    :param vec: vector to take the skew-symmetric matrix of
    :return skew: skew-symmetric matrix of vector vec
    """
    skew = np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])

    return skew


def euler_to_dcm(phi, theta, psi):
    """
    Return inertial-to-body transformation matrix from roll, pitch, and yaw Euler angles. Utilizes 3-2-1 rotation
    sequence. Reference [1].

    :param phi: roll (rad)
    :param theta: pitch (rad)
    :param psi: yaw (rad)
    :return dcm: inertial-to-body transformation matrix
    """
    dcm = rot_1(phi) @ rot_2(theta) @ rot_3(psi)

    return dcm
