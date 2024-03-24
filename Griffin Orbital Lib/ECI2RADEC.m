% Function to convert inertial position to right ascension and declination 
% 
% Author: Griffin Jourda 9/25/22
% 
%	Inputs 
%		r	:	inertial position vector 
%	Outputs 
%		ra	:	right ascension (rad) 
%		dec	:	declination (rad)

function [ra, dec] = ECI2RADEC(r)
	ra = atan2(r(2), r(1));
	dec = atan2(r(3), sqrt(r(1)^2 + r(2)^2));
end