% Function to convert a quaternion to a 3x3 rotation matrix 
% 
% Author: Griffin Jourda 11/16/22
% 
%	Inputs 
%		q	:	attitude quaternion [q1; q2; q3; qs]
% 
%	Outputs 
%		R	:	3x3 transformation matrix 
function [R] = q2R(q) 
	qs = q(4); 
	q = q(1:3);
	
	qx = [0, -q(3), q(2); 
		  q(3), 0, -q(1);
		  -q(2), q(1), 0];

	R = (qs^2 - q'*q)*eye(3) - 2*qs*qx + 2*(q*q');
end
