% Function to solve Lambert's problem with universal variables
% 
% Author: Griffin Jourda 10/12/22 (original) 9/30/24 (new version)
% 
% Inputs 
%	r1		:	initial position vector 
%	r2		:	final position vector 
%	tof		:	time of flight 
%	mu		:	central body gravitational parameter 
%	long_way:	long way (1) or short way (0) transfer
%	alpha	:	Newton-Raphson step size multiplier
% 
% Outputs 
%	v1	:	initial velocity 
%	v2	:	final velocity 
function [v1, v2] = lambert_uv(r1, r2, tof, mu, long_way, alpha)
	% Initial and final position mangitudes
	r1n = norm(r1); 
	r2n = norm(r2);
	
	% Change in true anomaly and direction of motion
	delta_f = acos(dot(r1, r2)/(r1n*r2n));
	DM = 1;
	if long_way
		delta_f = 2*pi - delta_f;
		DM = -1;
	end
	
	% Constant parameters
	A = DM*sqrt(r1n*r2n*(1 + cos(delta_f)));
	
	% Iterate over universal variable z
	err = 100;
	z = 0;
	while abs(err) > 1e-6
		% Universal variable computations
		[c, s] = stumpff(z);
		y = r1n + r2n - A*(1 - z*s)/sqrt(c);
		x = sqrt(y/c);
		
		% Time of flight
		t = (x^3*s + A*sqrt(y))/sqrt(mu);

		% Iteration
		[cprime, sprime] = stumpff_derivatives(z);
		dtdz = (x^3)*(sprime - 3*s*cprime/(2*c)) + (A/8)*(3*s*sqrt(y)/c + A/x);
		dtdz = dtdz/sqrt(mu);
		err = tof - t;
		z = z + alpha*err/dtdz;
	end
	
	% F and G functions
	F = 1 - y/r1n; 
	G = A*sqrt(y/mu);
	Gdot = 1 - y/r2n;

	% Velocities 
	v1 = (r2 - F*r1)/G; 
	v2 = (Gdot*r2 - r1)/G;
end