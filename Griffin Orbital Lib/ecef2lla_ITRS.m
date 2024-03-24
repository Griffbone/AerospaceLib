% Function to calculate latitude, longitude, and alittude from an ITRS
% position vector using the ITRS model. 
%
% Author: Griffin Jourda 3/18/2023
%
% Inputs
%	r	:	ITRF position vector (m) 
%
% Outputs
%	lat	:	latitude (deg) 
%	lon	:	longitude (deg)
%	alt	:	geodetic altitude (m)

function [lat, lon, alt] = ECEF2LLA_ITRS(r)
	% Constants for Earth (JGM-3 model reported in Vallado)
	Re = 6378136.3;
	ee = 0.081819221456;
		
	% Loop to calculate latitude
	rds = sqrt(r(1)^2 + r(2)^2); 
	d = atan2(r(3), rds); 
	lat0 = d;
	rd = rds;	
	err = 1;
	
	while err > 1e-9
		C = Re/sqrt(1 - ee^2*sin(lat0)^2);
		lat = atan2(r(3) + C*ee^2*sin(lat0), rd);
		err = abs(lat0 - lat); 
		lat0 = lat;
	end
	
	alt = rd/cos(lat) - C;
	lat = rad2deg(lat);
	lon = atan2d(r(2), r(1)); 
end