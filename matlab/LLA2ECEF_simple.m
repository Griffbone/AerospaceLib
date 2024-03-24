% Function to determine ECEF position vector from geocentric latitude,
% longitude and altitude. Assumes spherical Earth (R = 6378km); 
%
% Author: Griffin Jourda 10/04/22
% 
%	Inputs 
%		lat		:	geocentric latitude (deg)
%		lon		:	geocentric longitude (deg) 
%		alt		:	geocentric altitude (m) 
% 
%	Outputs 
%		r_ecef	:	ECEF position vector (m)

function [r_ecef] = LLA2ECEF_simple(lat, lon, alt)
	lat = deg2rad(lat); 
	lon = deg2rad(lon); 

	r = 6378e3 + alt;
	r_ecef = r*[cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];
end