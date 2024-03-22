import numpy as np
from rotationLib import rot_1, rot_2, rot_3


def lla_to_ecef(lat, lon, alt):
    """
    Return Earth-centered Earth-fixed (ECEF) coordinates from latitude, longitude, and altitude values.
    Utilizes WGS-84 values.

    :param lat: latitude (rad)
    :param lon: longitude (rad)
    :param alt: altitude (m)
    :return r_ecef: ECEF position vector (m)
    """
    re = 6378137
    ee = 0.081819221456
    ee2 = ee**2

    slat = np.sin(lat)
    clat = np.cos(lat)
    slon = np.sin(lon)
    clon = np.cos(lon)

    c = re/np.sqrt(1 - ee2*slat**2)
    s = re*(1 - ee2)/np.sqrt(1 - ee2*slat**2)

    x = (c + alt)*clat*clon
    y = (c + alt)*clat*slon
    z = (s + alt)*slat

    return np.array([x, y, z])
