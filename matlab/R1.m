% Create a counterclockwise rotation matrix about the X-axis
%
% Author: Griffin Jourda 2/11/2022
%
%	Inputs
%		theta	: rotation angle (rad)
%
%	Ouputs
%		rot		: rotation matrix

function [rot] = R1(theta)
	c = cos(theta); 
	s = sin(theta); 
	
	rot =	[1 0 0 
			0 c s
			0 -s c];
end