% Function to generate orbital elements of a sun-synchronous orbit about
% Earth. 
% 
% Author: Griffin Jourda 3/2/2023
%
% Inputs 
%	alt	:	orbital altitude (m) 
%	jd	:	Julian date at orbit insertion 
%	ltan:	local time of ascending node (hrs) 
% 
% Outputs 
%	r	:	ECI position vector (m) 
%	v	:	ECI velocity vector (m/s)

function [r, v] = generateSSO(alt, jd, ltan)
	re = 6378137;
	mu = 398600435436000;
	J2 = 1.08263e-3;

	a = re + alt;
	T = 2*pi*sqrt(a^3/mu);

	% Determine inclination
	rho = 1.99096871e-7;
	cf = @(i) ((-3*pi*J2*re^2*cos(i))/(a^2*T) - rho)*1e6;
	options = optimset('Display','off');
	i = fsolve(cf, 95*pi/180, options);

	% raan
	r_sun = sunVector(jd);
	raan = mod(atan2(r_sun(2), r_sun(1)) + ltan*2*pi/24, 2*pi);

	[r, v] = COE2RV(a, 0, i, raan, 0, 0, mu);
end