% Function to compute the equations of motion for translational and
% rotational state of a satellite, with control input and J2 perturbations.
% 
% Author: Griffin Jourda 9/25/22
%
%	Inputs
%		t			:	time (s)
%		y			:	state vector [r; v; q; w] (m, m/s, rad/s)
%		I			:	moment of inertia tensor (kg*m^2)
%		jd0			:	Julian date at simulation epoch (days)
%		controlfun	:	attitude controller function - must take time and
%						state as inputs and return a column vector of
%						moments in N*m about each body axis: 
%						M = controlfun(t, y)
%	Outputs 
%		ydot		:	state vector derivative [v; a; qdot; wdot] 
%						(m/s, m/s^2, rad/s^2)

function [ydot] = fullStateEOM(t, y, I, jd0, controlfun)
	mu = 398600435507000;
	J2 = 0.00108262672;
	Re = 6378135;

	% Unpack state vectors
	r = y(1:3); 
	v = y(4:6);
	q = y(7:10);
    q = q/norm(q);
	w = y(11:13);

	% Time conversions 
	jd = jd0 + t/86400;
	
	% J2 Acceleration
	rm = norm(r); 
	xi = r(1); 
	yi = r(2); 
	zi = r(3); 

	a_j2 = -(3/2)*J2*(mu/rm^2)*(Re/rm)^2;
	a_j2 = a_j2*[(xi/rm)*(1 - 5*(zi/rm)^2), (yi/rm)*(1 - 5*(zi/rm)^2), (zi/rm)*(3 - 5*(zi/rm)^2)];

	% Total acceleration 
	a = -(mu/rm^3)*r + a_j2';

	% Propagate rotational rate
	M = controlfun(t, y);
	wdot = I\(M - cross(w, I*w));
	qdot = (1/2)*qMult([w; 0], q); 
	
	% Compile state derivative 
	ydot = [v; a; qdot; wdot];
end