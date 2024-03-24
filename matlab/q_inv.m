% Function to perform a quaternion inverse 
% 
% Author: Griffin Jourda 11/16/22
% 
%	Inputs 
%		q		:	quaternion [q1; q2; q3; qs]
% 
%	Outputs
%		qinv	:	inverted quaternion 
function [qinv] = q_inv(q) 
	qinv = q;
	qinv(1:3) = -q(1:3); 
	qinv(4) = q(4); 
end
