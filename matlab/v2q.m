% Function to convert an angle vector to a quaternion 
% 
% Author: Griffin Jourda 11/15/22
% 
%	Inputs 
%		phi		:	3x1 angle vector (rad) 
% 
%	Outputs 
%		q		:	quaternion [q1; q2; q3; qs]
function [q] = v2q(phi) 
	theta = norm(phi); 
	e = phi/norm(phi); 
	q = [sin(theta/2)*e; cos(theta/2)];
end
