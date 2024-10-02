% Function to evaluate the derivative of the stumpff functions for
% universal variable z. For near-zero z values a series expression
% (truncated at four terms) is used. 
% 
% Author: Griffin Jourda 9/30/24
% 
% Inputs
%	z	:	universal variable z
% Outputs 
%	cprime	:	c' value 
%	sprime	:	s' value
function [cprime, sprime] = stumpff_derivatives(z) 
	[c, s] = stumpff(z);

	if abs(z) > 1e-6
		cprime = (1/(2*z))*(1 - z*s - 2*c);
		sprime = (1/(2*z))*(c - 3*s);
	else
		cprime = 0; 
		sprime = 0;
		for k = 1:4
			cprime = cprime + (((k+1)*z^k)/factorial(4 + 2*k))*(-1)^(k+1);
			sprime = sprime + (((k+1)*z^k)/factorial(5 + 2*k))*(-1)^(k+1);
		end
	end
end