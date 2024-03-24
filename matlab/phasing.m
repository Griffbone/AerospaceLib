% Function to calculate parameters for a phasing orbit from initial orbit,
% phase angle, and number of revolutions for phasing. Assumes circular
% initial orbit, elliptical phasing orbit, and two equal tangential burns. 
% 
% Author: Griffin Jourda 9/30/22
%
%	Inputs
%		a		:	initial semimajor axis
%		theta	:	angle to phase - positive eastward (rad)
%		nrev	:	number of revolutions to complete phasing
%		mu		:	gravitational parameter 
%
%	Outputs
%		a2		:	phasing orbit semimajor axis 
%		e		:	phasing orbit eccentricity 
%		T		:	total time to complete phasing
%		dv		:	total delta-V

function [a2, e, T, dv] = phasing(a, theta, nrev, mu)
	% Initial orbit
	n1 = sqrt(mu/a^3);
	T1 = 2*pi/n1;

	% Phasing orbit
	T2 = T1 - theta/(n1*nrev);
	a2 = (mu*(T2/(2*pi))^2)^(1/3);
	T = T2*nrev;

	% Delta-V
	v1 = sqrt(mu/a);
	v2 = sqrt(mu*(2/a - 1/a2)); 
	dv = 2*abs(v1 - v2); 

	h = v2*a;
	eps = (v2^2)/2 - mu/a;

	e = sqrt(1 + (2*eps*h^2)/mu^2);
end