% Function to signal an event to MATLAB ode solvers when a SOI is crossed 
% 
% Author: Griffin Jourda 10/13/22
% 
%	Inputs 
%		t	:	integration time (s) 
%		y	:	current state in the heliocentric J2000 ecliptic frame (km)
%		j0	:	Julian date at beginning of integration 
%		soi :	target sphere of influence radius 
%		pn	:	planet number (1 = Mercury, 2 = Venus, ... 8 = Neptune)
%
%	Outputs 
%		d			:	distance from the SOI 
%		isterminal	:	1 (stop integration when d == 0)
%		direction	:	0 (event can be approached on ascending or
%						descending d) 

function [d, isterminal, direction] = SOIEvent(t, y, j0, pn)	
	jd = j0 + t/86400;
	[rp, ~, mu] = planetState(jd, pn); 
	
	soi = norm(rp)*(mu/132712440018)^(2/5);
	
	r = y(1:3); 	
	d = norm(r' - rp) - soi;
	isterminal = 1; 
	direction = 0;
end