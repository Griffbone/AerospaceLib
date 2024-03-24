% Function to propagate an orbit based on a perturbed Kepler's problem.
% Function follows algorithm 65 in Vallado. 
%
% Author: Griffin Jourda 10/01/22
% 
%	Inputs 
%		r	:	initial position vector (m) 
%		v	:	initial velocity vector (m/s)
%		t	:	time span (s)
%		nd	:	mean motion derivative (rad/s^2)
%		ndd :	mean motion second derivative (rad/s^3)
%		mu	:	gravitational parameter (m^3/s^2)
% 
%	Outputs
%		t			:	time span (s)
%		y_out		:	array of satellite states (m, m/s)
%		elements	:	array of orbital elements [a, e, i, raan, w, ta] (m, rad)

function [t, y_out, elements] = keplerPropPerturbed(r, v, t, nd, ndd, mu)
	% Constants
	Re = 6378e3; 
	J2 = 0.0010826269; 
	
	% Initial elements
	[a0, e0, i0, raan0, w0, ta0] = RV2COE(r, v, mu);
	
	E0 = atan2(sin(ta0)*sqrt(1 - e0^2), e0 + cos(ta0));
	M0 = E0 - e0*sin(E0);
	n0 = sqrt(mu/a0^3);
	p0 = a0*(1 - e0^2);
	
	% Loop through time span 
	y_out = [r, v];
	elements = [a0, e0, i0, raan0, w0, ta0];
	
	for ti = t(2:end)		
		% Update orbital elements from perturbation 
		a = a0 - ((2*a0)/(3*n0))*nd*ti;
		e = e0 - (2*(1 - e0)/(3*n0))*nd*ti;
		raan = raan0 - ((3*n0*Re^2*J2)/(2*p0^2))*cos(i0)*ti;
		w = w0 + ((3*n0*Re^2*J2)/(4*p0^2))*(4 - 5*sin(i0)^2)*ti;
		
		% Update mean anomaly 
		M = M0 + n0*ti + (nd/2)*ti^2 + (ndd/6)*ti^3;
		
		% Solve Kepler's equation 
		[~, ta] = keplerE(M, e);
		[r, v] = COE2RV(a, e, i0, raan, w, ta, mu); 
		
		% Log data
		y_out = [y_out; r, v];
		elements = [elements; a, e, i0, raan, w, ta];
	end
end