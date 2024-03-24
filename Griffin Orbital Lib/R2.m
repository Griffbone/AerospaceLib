% Create a counterclockwise rotation matrix about the Y-axis
%
% Author: Griffin Jourda 2/11/2022
%
%	Inputs
%		theta	: rotation angle (rad)
%
%	Ouputs
%		rot		: rotation matrix

function [rot] = R2(theta)
	c = cos(theta); 
	s = sin(theta); 
	
	rot =	[c 0 -s;
			0 1 0;
			s 0 c];
end