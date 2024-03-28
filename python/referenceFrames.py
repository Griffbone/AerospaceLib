## @package referenceFrames
#  Library for common coordinate frame and reference system transformations used in aerospace engineering
#  
#  @author Griffin Jourda
import numpy as np
from mathematics import rot_1, rot_2, rot_3

## eci_to_ecef_gmst
#  @brief convert from ECI to ECEF coordinates with GMST transformation
#
#  This function computes ECEF position from ECI position using a simple conversion based on Greenwich Mean 
#  Sidereal Time (GMST)
#  @param r_eci ECI position vector (m) 
#  @param jd_ut1 UT1 Julian date (days)     # TODO: check which time system this should be in 
#  @return r_ecef ECEF position vector (m)
def eci_to_ecef_gmst(r_eci: np.array, jd_ut1: float): # TODO: check type hints here 
    # TODO: complete this function
    pass

## eci_to_ecef_iau
#  @brief convert from ECI to ECEF with full IAU transformation 
# 
#  This function computes ECEF position expressed in the ITRF frame from ECI position expressed in the ICRF frame.
#  The function utilizes the series method outlined in TODO: source/method year (IAU 1980?) and optionally takes 
#  pole position arugments for higher accuracy. 
#  @param r_eci ECI position vector in ICRF frame (m) 
#  @param jd_ut1 UT1 Julian date (days) TODO: check which time system this should be in
#  @param xp pole x position (m) 
#  @param yp pole y position (m)
def eci_to_ecef_iau(r_eci: np.array, jd_ut1: float, xp: float, yp: float):
    # TODO: write the rest of this function
    pass

## lla_to_ecef
# @brief convert from lat/lon to ECEF coordinates 
# 
# This function computes ECEF position from latitude and longitude expressed in the WGS-84 geodetic model. 
# @param lat geodetic latitude (rad) 
# @param lon geodetic longitude (rad) 
# @param alt altitude above WGS-84 ellipsoid (m)
# @return r_ecef ECEF position vector (m)
# @note source for WGS-84 constants: https://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
# @note discussion on WGS-84 and ITRS: https://geodesy.noaa.gov/CORS/Articles/Reference-Systems-Part-3.pdf
def lla_to_ecef(lat: float, lon: float, alt: float) -> np.array: # TODO: check type hints here
    a = 6378137.0
    f = 1/298.257223563
    ee = 2*f - f**2

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)

    c = a/np.sqrt(1 - ee*sin_lat**2)
    s = a*(1 - ee)/np.sqrt(1 - ee*sin_lat**2)

    x = (c + alt)*cos_lat*cos_lon
    y = (c + alt)*cos_lat*sin_lon
    z = (s + alt)*sin_lat

    return np.array([x, y, z])

## ecef_to_enu
# @brief convert from ECEF to ENU coordinates
# 
# This function computes position in an ENU frame centered around a ground station at a given latitude, longitude,
# and altitude
# @param r_ecef ECEF position vector (m)
# @param lat_gs geodetic ground station latitude (rad)
# @param lon_gs geodetic ground station longitude (rad)
# @param alt_gs ground station altitude above WGS-84 ellipsoid (m)
def ecef_to_enu(r_ecef: np.array, lat_gs: float, lon_gs: float, alt_gs: float) -> np.array: # TODO: check type hints here
    # ground station ECEF position
    r_gs_ecef = lla_to_ecef(lat_gs, lon_gs, alt_gs);

    # line of sight vector
    s = r_ecef - r_gs_ecef

    # rotate to ENU
    r_enu = rot_1(np.pi/2 - lat_gs) @ rot_3(np.pi/2 + lon_gs) @ s
    
    return r_enu

## enu_to_azel
# @brief convert from ENU coordinates to azimuth and elevation 
# 
# This function computes azimuth and elevation from the origin of an ENU frame to a position 
# vector expressed within that frame
# @param 
def enu_to_azel(r_enu: np.array) -> tuple[float, float]: # TODO: check type hints on np array
    l = np.sqrt(r_enu[0]**2 + r_enu[1]**2)
    el = np.arctan2(r_enu[3], l)
    az = 0  # TODO: put in actual computation here

    return az, el