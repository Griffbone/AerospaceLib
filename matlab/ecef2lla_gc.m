% Function to calculate geocentric latitude, longitude, and altitude from
% an ECEF position vector. Assumes spherical Earth with r = 6378 km. 
%
% Author: Griffin Jourda 3/18/22
% 
%	Inputs
%		r	:	ECEF position vector (m) 
% 
%	Outputs 
%		lat	:	geocentric or geodetic latitude (rad)
%		lon	:	geocentric or geodetic longitude (rad)
%		alt :	geocentric or geodetic altitude (m)

function [lat, lon, alt] = ecef2lla_gc(r)
	alt = sqrt(r(:, 1).^2 + r(:, 2).^2 + r(:, 3).^2) - 6378e3;
	lon = atan2(r(:, 2), r(:, 1));
	lat = atan2(r(:, 3), sqrt(r(:, 1).^2 + r(:, 2).^2));
end