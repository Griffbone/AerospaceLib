% Calculate Greenwich mean sidereal time at given Julian date
%
% Author: Griffin Jourda 9/27/2022
%
%	Inputs
%		jd		:	Julian date 
%
%	Ouputs
%		theta_g	:	Greenwich mean sidereal time (rad)

function [theta_g] = jd2gmst(jd)
	% Centuries since J2000 
	T = (jd - 2451545.0)/36525;

	% Greenwich mean sidereal time (s)
	theta_g = 67310.54841 + (876600*60^2 + 8640184.812866)*T + 0.09304*T^2 - 6.2e-6*T^3;

	% Wrap to 360 and convert to radians 
	theta_g = deg2rad(mod(theta_g/240, 360));
end