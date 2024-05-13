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
%		r_ecef		:	position vector in ECEF frame
%		r_eci_ecef	:	DCM from ECI to ECEF frame

function [r_ecef, r_eci_ecef] = eci2ecef_gmst(r, jd)
	era = jd2gmst(jd); 
	r_eci_ecef = R3(era);
	r_ecef = r_eci_ecef*r;
end