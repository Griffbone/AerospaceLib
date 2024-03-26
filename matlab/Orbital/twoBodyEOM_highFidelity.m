% Function to calculate the derivative of a positon, velocity state vector
% for an Earth satellite, taking into account pertubations from
% nonspherical gravity, lunisolar perturbations, solar radiation pressure,
% and atmospheric drag. Utilizes low-precision Sun and Moon ephemerides and
% the Harris-Priester drag model.
% 
% Author: Griffin Jourda 4/5/2023
% 
% Inputs
%	t			:	current time (s)
%	y			:	current state [m, m, m, m/s, m/s, m/s]
%	jd0			:	Julian date at start of simulation 
%	m			:	satellite mass (kg)
%	cda_drag	:	satellite drag area times drag coefficient
%	cra_srp		:	satellite srp coefficient times srp area 
%	EOPData		:	Earth orientation parameter data for interpolation for
%					ECI-->ECEF reduction [jd, xp, yp, UT1-UTC, UT1-TAI]
%	C			:	array of denormalized gravity model C coefficients 
%	S			:	array of denormalized gravity model S coefficients 
%	mu			:	Earth gravitational parameter for utilized gravity
%					model
%	rp			:	Earth radius for utilized gravity model
% Outputs
%	ydot		:	state derivative [m/s, m/s, m/s, m/s^2, m/s^2, m/s^2]

function [ydot] = twoBodyEOM_highFidelity(t, y, jd0, m, cda_drag, cra_srp, EOPData, C, S, mu, rp)
	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);

	% Time and EOP parameters
	jd_utc = jd0 + t/86400;
	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
	
	% Coordinate transformations
	R_icrf_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);
	r_ecef = R_icrf_itrf*r_eci;
	
	% Earth gravity 
	[gx, gy, gz] = sphericalHarmonicGravity(r_ecef, C, S, rp, mu);
	a_earth = (R_icrf_itrf')*[gx; gy; gz];
	
	% Moon/Sun gravity
	s_sun = sunVector(jd_utc);
	smr_sun = s_sun - r_eci;
	s_moon = moonVector(jd_utc);
	smr_moon = s_moon - r_eci;

	a_sun = 1.32712440041939e20*((smr_sun)/norm(smr_sun)^3 - s_sun/norm(s_sun)^3);
	a_moon = 4902800066000*((smr_moon)/norm(smr_moon)^3 - s_moon/norm(s_moon)^3);
	
	% Solar radiation pressure
	if ~isEclipsed(r_eci, s_sun)
		r_sun = s_sun;
		f_srp = -(-4.56e-6*cra_srp)*(r_sun/norm(r_sun)^3)*2.23795229152812e22;
	else 
		f_srp = 0;
	end
	
	% Drag
	h = cross(r_eci, v_eci);
	inc = acos(h(3)/norm(h));
	[~, ~, alt] = ecef2lla_ITRS(r_ecef);
	rho = harrisPriester(alt, r_eci, s_sun, inc);
	v_rel = v_eci - cross([0; 0; 2*pi/86164], r_ecef);
	f_d = -(1/2)*cda_drag*rho*norm(v_rel)*v_rel;

	% Total state derivative 
	a = a_earth + a_sun + a_moon + (f_d + f_srp)*(1/m);
	ydot = [v_eci; a];
end