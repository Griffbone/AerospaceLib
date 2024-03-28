import numpy as np
from python.mathematics import euler_to_dcm, skew_sym


def eom(_, y, f_b, m_b, m, I):
    """
    Function to return derivative of aircraft motion state in body-frame coordinates using roll-pitch-yaw Euler angles.

    :param y: aircraft state [xyz, uvw, rpy, pqr] (m, m/s, rad, rad/s)
    :param f_b: body-frame forces on the aircraft (N)
    :param m_b: body-frame moments on the aircraft (N*m)
    :param m: mass (kg)
    :param I: aircraft inertia tensor (kg*m**2)
    :return y_dot: state derivative [xyz_dot, uvw_dot, rpy_dot, pqr_dot] (m/s, m/s**2, rad/s, rad/s**2)
    """
    v_b = y[3:6]
    rpy = y[6:9]
    w_b = y[9:12]

    # Translational kinematics
    dcm_i_b = euler_to_dcm(rpy[0], rpy[1], rpy[2])
    dcm_b_i = dcm_i_b.T
    r_i_dot = dcm_b_i @ v_b

    # Euler angle kinematics
    s_phi = np.sin(rpy[0])
    c_phi = np.cos(rpy[0])
    c_theta = np.cos(rpy[1])
    t_theta = np.tan(rpy[1])

    rpy_dot = np.array([[1, s_phi*t_theta, c_phi*t_theta],
                        [0, c_phi, -s_phi],
                        [0, s_phi/c_theta, c_phi/c_theta]]) @ w_b

    # Translational dynamics
    w_b_skew = skew_sym(w_b)
    v_b_dot = (f_b/m) - w_b_skew @ v_b

    # Rotational dynamics
    w_b_dot = np.linalg.inv(I) @ (m_b - w_b_skew @ I @ w_b)

    # Compile state derivative
    y_dot = np.concatenate((r_i_dot, v_b_dot, rpy_dot, w_b_dot))

    return y_dot
