% Function to calculate geocentric position of the Sun from Julian date.
% Function follows "simple" algorithm laid out in Meeus (1998) ch. 25. 
% 
% Author: Griffin Jourda 10/06/22
% 
%	Inputs
%		jd	:	Julian date		
%		
%	Outputs
%		ra	:	apparent right ascension of the sun (deg) 
%		dec :	apparent declination of the sun (deg) 

function [ra, dec] = sunRaDec(jd) 
	clc

	T = (jd - 2451545.0)/36525;

	L0 = 280.46646 + 36000.76983*T + 0.0003032*T^2;
	M = 357.52911 + 35999.05029*T - 0.0001537*T^2;
	e = 0.016708634 - 0.000042037*T - 0.0000001267*T^2;

	C = (1.914602 - 0.004817*T - 0.000014*T^2)*sind(M) + ...
		(0.019993 - 0.000101*T)*sind(2*M) + 0.000289*sind(3*M);
	
	% True longitude - mean equinox of date 
	trueLon = L0 + C; 
	trueLon = mod(trueLon, 360);

	ta = M + C; 
	R = (1.000001018*(1 - e^2))/(1 + e*cosd(ta));
	Omega = 125.04 - 1934.136*T; 

	% Apparent longitude - true equinox of date 
	lon = trueLon - 0.00569 - 0.00478*sind(Omega);

	% Obliquity of ecpliptic, Ra, and Dec
	eps = (23 + 26/60 + 21.448/3600) - (46.8150/3600)*T - (0.00059/3600)*T^2 + (0.001813/3600)*T^2; 
	eps = eps + 0.00256*cosd(Omega); 
	ra = atan2d(cosd(eps)*sind(lon), cosd(lon));
	dec = asind(sind(eps)*sind(lon));
	ra = mod(ra, 360);
end