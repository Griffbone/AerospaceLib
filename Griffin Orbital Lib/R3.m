% Create a counterclockwise rotation matrix about the Z-axis
%
% Author: Griffin Jourda 2/11/2022
%
%	Inputs
%		theta	: rotation angle (rad)
%
%	Ouputs
%		rot		: rotation matrix

function [rot] = R3(theta)
	c = cos(theta); 
	s = sin(theta); 
	
	rot =	[c s 0;
			-s c 0;
			0 0 1];
end