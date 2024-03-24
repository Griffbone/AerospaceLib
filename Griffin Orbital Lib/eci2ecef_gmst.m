% Function to convert from inertial ECI frame to rotating ECEF frame. This
% function solely uses a z-rotation by the Earth rotation angle and does
% not account for procession, nutation, or polar motion.
%
% Author: Griffin Jourda 9/25/22
% 
%	Inputs 
%		r		:	position vector in ECI frame (m) 
%		jd		:	UT Julian date 
%
%	Outputs 
%		r_ecef	:	position vector in ECEF frame

function [r_ecef] = ECI2ECEF_simple(r, jd)
	era = jd2gmst(jd); 
	r_ecef = R3(era)*r;
end