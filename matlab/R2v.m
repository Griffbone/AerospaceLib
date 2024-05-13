% Function to convert a rotation matrix into an angle vecotr 
% 
% Author: Griffin Jourda 11/30/22
% 
%	Inputs 
%		R	:	3x3 rotation matrix 
%	
%	Outputs 
%		v	:	angle vector (rad)
function [v] = R2v(R)
	[v, d] = eig(R); 
	d = sum(real(d));
	i = abs(d-1) == min(abs(d-1));
	v = v(:, i); 

	theta = acos((trace(R) - 1)/2);
	v = v*theta; 
end 