% Function to return derivative of the state vector for a satellite moving
% around Earth under two-body motion 
% 
% Author: Griffin Jourda 9/25/22
%
%	Inputs
%		t		:	time (s)
%		y		:	state vector [rx; ry; rz; vx; vy; vz;] (m, m/s)
%	Outputs 
%		ydot	:	state vector derivative [vx; vy; vz; ax; ay; az] (m, m/s^2)

function [ydot] = twoBodyEOM(~, y, mu)
	r = y(1:3); 
	v = y(4:6);

	a = -(mu/norm(r)^3)*r;
	ydot = [v; a];
end