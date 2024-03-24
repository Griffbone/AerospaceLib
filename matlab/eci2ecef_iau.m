% Function to calculate the rotation matrix from the International
% Celestial Reference Frame to the International Terrestrial Reference
% Frame. Utilizes IAU-1980 nutation model. Follows algorithms in Montebruck
% and Gill, Chapter 5. 
%
% Inputs 
%	jd_tt		:	Terrestrial time Julian date
%	jd_ut1		:	UT1 Julian date
%	xp			:	X pole position (arcseconds)
%	yp			:	Y pole position (arcseconds)
%
% Outputs 
%	R_icrf_itrf		:	Rotation matrix from ICRF to ITRF 
%	R_icrf_itrf_p	:	Derivative of the rotation matrix

function [R_icrf_itrf, R_icrf_itrf_p] = eci2ecef_iau(jd_tt, jd_ut1, xp, yp)
	T_tt = (jd_tt - 2451545)/36525;
	T_tt2 = T_tt^2; 
	T_tt3 = T_tt^3; 
	D2R = pi/180;
	
	% Precession
	zeta = (2306.2181/3600)*T_tt + (0.30188/3600)*T_tt2 + (0.017998/3600)*T_tt3;
	nu = (2004.3109/3600)*T_tt - (0.42665/3600)*T_tt2 - (0.041833/3600)*T_tt3;
	z = zeta + (0.79280/3600)*T_tt2 + (0.000205/3600)*T_tt3;
	P = R3(-z*D2R)*R2(nu*D2R)*R3(-zeta*D2R);
	
	% Nutation 
	[eps, dpsi, deps] = nutationIAU1980(jd_tt);
	N = R1((-eps - deps))*R3((-dpsi))*R1(eps);
	
	% Earth rotation 
	GMST = jd2gmst(jd_ut1);
	GAST = GMST + dpsi*cos(eps);
	Theta = R3(GAST);

	dThetadt = (7.2921158553e-5)*[0 1 0; -1 0 0; 0 0 0]*Theta;
	
	% Polar motion
	xp = (xp/3600)*D2R;
	yp = (yp/3600)*D2R;
	Pi = [1 0 xp; 0 1 -yp; -xp yp 1];
	
	% Combined transformation
	R_icrf_itrf = Pi*Theta*N*P;
	R_icrf_itrf_p = Pi*dThetadt*N*P;
end