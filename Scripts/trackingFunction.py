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
    el = np.arctan2(r_enu[2], l)
    az = np.arctan2(r_enu[0], r_enu[1]) % (2*np.pi)

    return az, el

## lla_to_azel 
# @brief convert a rocket position to antenna azimuth and elevation 
# 
# This functin receives a longitude and latitude representing the position 
# of a vehicle and returns an azimuth (clockwise from true north) and elevation 
# (up from horizontal) representing the vector from a given ground station to the target.
# @param lat geodetic latitude (rad) 
# @param lon geodetic longitude (rad) 
# @param alt geodetic altitude (m) 
# @param lat_gs ground station geodetic latitude (rad) 
# @param lon_gs ground station geodetic longitude (rad) 
# @param alt_gs ground station geodetic altitude (m)
def lla_to_azel(lat: float, lon: float, alt: float,
                lat_gs: float, lon_gs: float, alt_gs: float) -> tuple[float, float]:
    # Compute vehicle ECEF position
    r_ecef = lla_to_ecef(lat, lon, alt)

    # Convert vehicle position to ground station ENU frame 
    r_enu = ecef_to_enu(r_ecef, lat_gs, lon_gs, alt_gs)

    # Convert vehicle enu position to azimuth and elevation
    az, el = enu_to_azel(r_enu)

    return az, el

# Example use 
lat_gs = 38*np.pi/180
lon_gs = -100*np.pi/180
alt_gs = 100

lat_rocket = 38*np.pi/180
lon_rocket = -100*np.pi/180
alt_rocket = 10e3

az, el = lla_to_azel(lat_rocket, lon_rocket, alt_rocket, lat_gs, lon_gs, alt_gs)
print(az*180/np.pi)
print(el*180/np.pi)