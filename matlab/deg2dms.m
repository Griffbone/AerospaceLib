% Function to convert an angle in degrees to degrees, minutes, seconds
% 
% Author: Griffin Jourda 10/06/22
% 
%	Inputs
%		deg	:	angle in degrees 
%	
%	Outputs 
%		d	:	degrees
%		m	:	minute angle 
%		s	:	second angle 

function [d, m, s] = deg2dms(deg)
	sign1 = sign(deg); 
	deg = abs(deg); 

	d = floor(deg); 
	m = floor((deg - d)*60); 
	s = (deg - d - m/60)*3600;

	d = sign1*d; 
end