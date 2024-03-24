% Function to evaluate the Stumpff functions
% 
% Author: Griffin Jourda 10/12/22
% 
%	Inputs
%		z	:	z value
% 
%	Outputs 
%		c	:	c value 
%		s	:	s value

function [c, s] = stumpff(z)
	if z < 0 
		c = (1 - cosh(sqrt(-z)))/z; 
		s = (sinh(sqrt(-z)) - sqrt(-z))/sqrt((-z)^3);
	elseif z > 0
		c = (1 - cos(sqrt(z)))/z; 
		s = (sqrt(z) - sin(sqrt(z)))/sqrt(z^3);
	else 
		c = 1/2; 
		s = 1/6; 
	end
end