% Function to convert a position vector in the ECEF frame to a slant range
% vector in the topocentric ENU frame. Assumes spherical earth (R = 6378
% km) and geocentric coordinates.
% 
% Author: Griffin Jourda 10/04/22
% 
%	Inputs
%		r_ecef	:	position vector in the ECEF frame (m)
%		lat		:	latitude of the topocentric site (deg) 
%		lon		:	longitude of the topocentric site (deg)
% 
%	Outputs 
%		r_enu	:	position vector in the topocentric ENU frame 

function [r_enu] = ECEF2ENU_simple(r_ecef, lat, lon)
	lat = deg2rad(lat); 
	lon = deg2rad(lon);

	R = [-sin(lon), cos(lon), 0; 
		-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat); 
		cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];

	r_enu = R*r_ecef;
	r_enu = (r_enu - [0 0 6378e3]');
end