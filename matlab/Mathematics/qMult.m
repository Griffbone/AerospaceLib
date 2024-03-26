% Function to multiply two quaternions 
% 
% Author: Griffin Jourda 11/16/22
% 
%	Inputs
%		q1	:	first quaternion [q1; q2; q3; qs]
%		q2	:	second quaternion [q1; q2; q3; qs]
%
%	Outputs 
%		q_out	:	multiplied quaternion [q1; q2; q3; qs]
function [q_out] = qMult(q1, q2)
	qa = q1(1:3); 
	qas = q1(4); 
	qb = q2(1:3); 
	qbs = q2(4); 

	q_out = [qbs*qa + qas*qb - cross(qa, qb); 
			 qas*qbs - qa'*qb];
end