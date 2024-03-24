% Function to convert a rotation matrix into an angle vecotr 
% 
% Author: Griffin Jourda 11/30/22
% 
%	Inputs 
%		R	:	3x3 rotation matrix 
%	
%	Outputs 
%		ax	:	euler axis
%		ang	:	rotation angle (rad)
function [ax, ang] = R2v(R)
	[v, d] = eig(R); 
	d = sum(real(d));
	i = abs(d-1) == min(abs(d-1));
	ax = v(:, i); 

	ang = acos((trace(R) - 1)/2);
end 