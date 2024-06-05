% Function to convert a position vector in the ECEF frame to a slant range
% vector in the topocentric ENU frame given positionin ECEF and ground
% station geodetic coordinates. 
% 
% Author: Griffin Jourda 5/22/2024
% 
%	Inputs
%		r_ecef	    :	position vector in the ECEF frame (m)
%		lat_gs		:	latitude of the topocentric site (rad) 
%		lon_gs		:	longitude of the topocentric site (rad)
%       alt_gs      :   altitude of the topocentric site (m)
% 
%	Outputs 
%		r_enu	:	position vector in the topocentric ENU frame 

function [r_enu] = ECEF2ENU(r_ecef, lat_gs, lon_gs, alt_gs)
    r_gs = LLA2ECEF(lat_gs, lon_gs, alt_gs);

	R = [-sin(lon_gs), cos(lon_gs), 0; 
		-sin(lat_gs)*cos(lon_gs), -sin(lat_gs)*sin(lon_gs), cos(lat_gs); 
		cos(lat_gs)*cos(lon_gs), cos(lat_gs)*sin(lon_gs), sin(lat_gs)];
    
    s = r_ecef - r_gs'; 
    r_enu = R*s;
end