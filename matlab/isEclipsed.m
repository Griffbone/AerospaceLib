% Function to determine if a satellite in an Earth orbit is eclipsed by
% Earth 
% 
% Author: Griffin Jourda 4/7/2023
% 
% Inputs 
%	r			:	geocentric satellite position (m) 
%	s			:	geocentric sun position (m) 
% Outputs
%	eclipse		:	logical 1 if sun is eclipsed, otherwise logical 0

function [eclipse] = isEclipsed(r, s) 
	theta = mod(asin(6378000/norm(r)), 2*pi);
	psi = mod(acos(dot(r, s)/(norm(r)*norm(s))), 2*pi);

	if psi < theta 
		eclipse = true;
	else 
		eclipse = false;
	end
end