% Function to calculate the skew-symmetric (cross-product) matrix of a vector 
% 
% Author: Griffin Jourda 10/16/23 
% 
% Inputs
%	a	:	3x1 vector 
% Outputs 
%	askew	:	skew-symmetric matrix of a
function [askew] = skewSym(a)
	askew = [0, -a(3), a(2); 
		a(3), 0, -a(1); 
		-a(2), a(1), 0];
end