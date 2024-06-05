% Function to determine ECEF position vector from geodetic latitude,
% longitude, and altitude. Utilized JGM-3 Earth model. 
%
% Author: Griffin Jourda 10/01/22
% 
%	Inputs 
%		lat		:	geodetic latitude (rad)
%		lon		:	geodetic longitude (rad) 
%		alt		:	geodetic altitude (m) 
% 
%	Outputs 
%		r_ecef	:	ECEF position vector (m)

function [r_ecef] = LLA2ECEF(lat, lon, alt)
	% Constants for Earth (JGM-3 model reported in Vallado)
% 	Re = 6378136.3;
	Re = 6378137;
	ee = 0.081819221456;
	
	C = Re/sqrt(1 - ee^2*sin(lat)^2);
	S = (Re*(1 - ee^2))/sqrt(1 - ee^2*sind(lat)^2);
	
	x = (C + alt)*cos(lat)*cos(lon); 
	y = (C + alt)*cos(lat)*sin(lon); 
	z = (S + alt)*sin(lat); 
	
	r_ecef = [x, y, z];
end