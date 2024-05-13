% Function to return derivative of the state vector for a satellite moving
% around Earth under two-body motion with J2 perturbations 
% 
% Author: Griffin Jourda 9/25/22
%
%	Inputs
%		t		:	time (s)
%		y		:	state vector [rx; ry; rz; vx; vy; vz;] (m, m/s)
%	Outputs 
%		ydot	:	state vector derivative [vx; vy; vz; ax; ay; az] (m, m/s^2)

function [ydot] = twoBodyEOM_J2(~, y) 
	mu = 398600435507000;
	J2 = 0.00108262672;
	Re = 6378135;

	rv = y(1:3); 
	r = norm(rv); 
	v = y(4:6);
	
	% J2 Acceleration
	x = rv(1); 
	y = rv(2); 
	z = rv(3); 

	a_j2 = -(3/2)*J2*(mu/r^2)*(Re/r)^2;
	a_j2 = a_j2*[(x/r)*(1 - 5*(z/r)^2), (y/r)*(1 - 5*(z/r)^2), (z/r)*(3 - 5*(z/r)^2)];

	% Total acceleration 
	a = -(mu/r^3)*rv + a_j2';

	ydot = [v; a];
end