% Function to find the approximate location of the sun in the J2000
% Earth-centered inertial frame. Function follows process given in section
% 3.2.2 of Montenbruck and Gill. 
%
% Author: Griffin Jourda 3/7/2023
% 
%	Inputs 
%		jd		:	UT Julian date 
%
%	Outputs 
%		r_sun	:	sun vector in ECI J2000 frame (m)

function [r_sun] = sunVector(jd)
	T = (jd - 2451545)/36525;

	om = 282.9400;				% deg
	M = 357.5256 + 35999.049*T;	% deg 
	eps = 23.43929111;			% deg

	lo = om + M + (6892/3600)*sind(M) + (72/3600)*sind(2*M);	% deg
	ro = (149.619 - 2.499*cosd(M) - 0.021*cosd(2*M))*10^9;		% m
	
	r_sun = [ro*cosd(lo);
			ro*sind(lo)*cosd(-eps); 
			-ro*sind(lo)*sind(-eps)];
end