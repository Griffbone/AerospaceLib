clc; clear;
jd0 = greg2jd(2023, 1, 1, 0, 0, 0);
mjd0_gmat = jd0 - 2430000.0;
mu_egm2008 = 0.3986004415e15;

[r0, v0] = generateSSO(550e3, jd0, 18);
fprintf('Initial Pos: %.6f, %.6f, %.6f km\n', r0(1)/1000, r0(2)/1000, r0(3)/1000)
fprintf('Initial Vel: %.6f, %.6f, %.6f km/s\n', v0(1)/1000, v0(2)/1000, v0(3)/1000)
fprintf('Initial GMAT MJD: %.6f \n', mjd0_gmat)

[r0, v0] = COE2RV(42164*1000, 0.0001, 0.02*pi/180, 0, 0, 0, mu_egm2008);
fprintf('Initial Pos: %.6f, %.6f, %.6f km\n', r0(1)/1000, r0(2)/1000, r0(3)/1000)
fprintf('Initial Vel: %.6f, %.6f, %.6f km/s\n', v0(1)/1000, v0(2)/1000, v0(3)/1000)
fprintf('Initial GMAT MJD: %.6f \n', mjd0_gmat)

%% Earth as a point mass only
clc; clear;
[t_gmat, y_gmat, ~] = interpGMAT('geo_twoBody.e');

options = odeset('InitialStep', 60, 'RelTol', 10e-12, 'MaxStep', 2700);
[t, y] = ode89(@(t, y) twoBodyEOM_local(t, y), linspace(0, max(t_gmat), length(t_gmat)), 1000*y_gmat(1, :)', options);

hold on
plot(t, y_gmat*1000 - y);
plot(t, vecnorm(y_gmat*1000 - y, 2, 2))

%% Earth, Sun, Moon Point Masses
clc; clear;
[t_gmat, y_gmat, jd0] = interpGMAT('ephem/sso_twoBodyLunisolar.e');

options = odeset('InitialStep', 60, 'RelTol', 10e-12, 'MaxStep', 2700);
[t, y] = ode89(@(t, y) twoBodyEOM_lunisolar(t, y, jd0), linspace(0, max(t_gmat), length(t_gmat)), 1000*y_gmat(1, :)', options);

hold on
plot(t, y_gmat*1000 - y);
plot(t, vecnorm(y_gmat*1000 - y, 2, 2))

%% Earth, Sun, Moon Point Masses, SRP
clc; clear;
[t_gmat, y_gmat, jd0] = interpGMAT('ephem/geo_twoBodyLunisolarSRP.e');

options = odeset('InitialStep', 60, 'RelTol', 10e-12, 'MaxStep', 2700);

m = 1000; 
cra_srp = 10*1.3;
[t, y] = ode89(@(t, y) twoBodyEOM_lunisolarSRP(t, y, jd0, m, cra_srp), linspace(0, max(t_gmat), length(t_gmat)), 1000*y_gmat(1, :)', options);

hold on
plot(t, y_gmat(:, 1:3)*1000 - y(:, 1:3));
plot(t, vecnorm(y_gmat(:, 1:3)*1000 - y(:, 1:3), 2, 2))

%% Earth Sun Moon point masses, Drag
clc; clear;
[t_gmat, y_gmat, jd0] = interpGMAT('ephem/sso_twoBodyLunisolarDrag.e');
EOPData = loadEOPData('eopData.txt');

options = odeset('InitialStep', 60, 'RelTol', 10e-12, 'MaxStep', 2700);

m = 1000; 
cda_drag = 10*2.2;
[t, y] = ode89(@(t, y) twoBodyEOM_lunisolarDrag(t, y, jd0, m, cda_drag, EOPData), linspace(0, max(t_gmat), length(t_gmat)), 1000*y_gmat(1, :)', options);

hold on
plot(t, y_gmat(:, 1:3)*1000 - y(:, 1:3));
plot(t, vecnorm(y_gmat(:, 1:3)*1000 - y(:, 1:3), 2, 2))

%% Earth nonspherical gravity
clc; clear;
[t_gmat, y_gmat, jd0] = interpGMAT('ephem/sso_twoBodySHGrav.e');
EOPData = loadEOPData('eopData.txt');
% [C, S] = loadSHData('EGM2008.txt', 3);

options = odeset('InitialStep', 60, 'RelTol', 10e-12, 'MaxStep', 2700);

[t, y] = ode89(@(t, y) twoBodyEOM_shgrav(t, y, jd0, EOPData, 0, 0), linspace(0, max(t_gmat), length(t_gmat)), 1000*y_gmat(1, :)', options);

hold on
plot(t, y_gmat(:, 1:3)*1000 - y(:, 1:3));
plot(t, vecnorm(y_gmat(:, 1:3)*1000 - y(:, 1:3), 2, 2))

%% Restricted force model functions
% Two-body equations of motion with only Earth point mass 
function [ydot] = twoBodyEOM_local(~, y)
	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);
	
	% Earth gravity 
	a = -(0.3986004415e15/norm(r_eci)^3)*r_eci;
	
	% Total state derivative 
	ydot = [v_eci; a];
end

% Two-body equations of motion with lunisolar perturbations 
function [ydot] = twoBodyEOM_lunisolar(t, y, jd0)
	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);
	jd_utc = jd0 + t/86400;
	
	% Earth gravity 
	a_earth = -(0.3986004415e15/norm(r_eci)^3)*r_eci;
	
	% Moon/Sun gravity
	s_sun = sunVector(jd_utc);
	smr_sun = s_sun - r_eci;
	s_moon = moonVector(jd_utc);
	smr_moon = s_moon - r_eci;

	a_sun = 1.32712440041939e20*((smr_sun)/norm(smr_sun)^3 - s_sun/norm(s_sun)^3);
	a_moon = 4.9048695e12*((smr_moon)/norm(smr_moon)^3 - s_moon/norm(s_moon)^3);
	
	% Total state derivative 
	a = a_earth + a_sun + a_moon;
	ydot = [v_eci; a];
end

% Two-body equations of motion with lunisolar perturbations and SRP
function [ydot] = twoBodyEOM_lunisolarSRP(t, y, jd0, m, cra_srp)
	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);
	jd_utc = jd0 + t/86400;
	
	% Earth gravity 
	a_earth = -(0.3986004415e15/norm(r_eci)^3)*r_eci;
	
	% Moon/Sun gravity
	s_sun = sunVector(jd_utc);
	smr_sun = s_sun - r_eci;
	s_moon = moonVector(jd_utc);
	smr_moon = s_moon - r_eci;

	a_sun = 1.32712440041939e20*((smr_sun)/norm(smr_sun)^3 - s_sun/norm(s_sun)^3);
	a_moon = 4.9048695e12*((smr_moon)/norm(smr_moon)^3 - s_moon/norm(s_moon)^3);

	% SRP
	if ~eclipseFinding(r_eci, s_sun)
		f_srp = -4.56e-6*cra_srp*(s_sun/norm(s_sun)^3)*2.23795229152812e22;
	else 
		f_srp = 0;
	end

	% Total state derivative 
	a = a_earth + a_sun + a_moon + (1/m)*(f_srp);
	ydot = [v_eci; a];
end

% Two-body equations of motion with lunisolar perturbations and drag
function [ydot] = twoBodyEOM_lunisolarDrag(t, y, jd0, m, cda_drag, EOPData)
	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);
	jd_utc = jd0 + t/86400;

	% Coordinate transformations
	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
	R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);
	r_ecef = R_icrs_itrf*r_eci;
	[~, ~, h_ellipsoid] = ecef2lla_ITRS(r_ecef);
	
	% Earth gravity 
	a_earth = -(0.3986004415e15/norm(r_eci)^3)*r_eci;
	
	% Moon/Sun gravity
	s_sun = sunVector(jd_utc);
	smr_sun = s_sun - r_eci;
	s_moon = moonVector(jd_utc);
	smr_moon = s_moon - r_eci;

	a_sun = 1.32712440041939e20*((smr_sun)/norm(smr_sun)^3 - s_sun/norm(s_sun)^3);
	a_moon = 4.9048695e12*((smr_moon)/norm(smr_moon)^3 - s_moon/norm(s_moon)^3);

	% Drag
	h = cross(r_eci, v_eci);
	inc = acos(h(3)/norm(h));
	rho = harrisPriester(h_ellipsoid, r_eci, s_sun, inc);
	v_rel = v_eci - cross([0; 0; 2*pi/86164], r_ecef);
	f_d = (-1/2)*cda_drag*rho*v_rel*norm(v_rel);

	% Total state derivative 
	a = a_earth + a_sun + a_moon + (1/m)*(f_d);
	ydot = [v_eci; a];
end

% Two-body equations of motion with nonspherical gravity
function [ydot] = twoBodyEOM_shgrav(t, y, jd0, EOPData, C, S)
	disp(t/86400);

	% Unpack state vector
	r_eci = y(1:3);
	v_eci = y(4:6);
	jd_utc = jd0 + t/86400;

	% Coordinate transformations
	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
	R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);
	r_ecef = R_icrs_itrf*r_eci;
	
	% Earth gravity 
% 	[x, y, z] = sphericalHarmonicGravity(r_ecef, C, S, 6378135, 0.3986004415e15);
	[x, y, z] = gravitysphericalharmonic(r_ecef', 'EGM96', 2);
	a_earth = R_icrs_itrf'*[x; y; z];
	
	% Total state derivative 
	a = a_earth; 
	ydot = [v_eci; a];
end







function [shade] = shadow(r, s)
	alf_umb = 0.264121687*pi/180; 
	alf_pen = 0.269007205*pi/180;
	shade = 0;

	if dot(r, s) < 0
		zeta = acosd(dot(-s, r)/(norm(-s)*norm(r)));

		horiz = norm(r)*cos(zeta);
		vert = norm(r)*sin(zeta); 
		x = 6378000/sin(alf_pen);

		penvert = tan(alf_pen)*(x + horiz);

		if vert < penvert 
			shade = 1;
			y = 6378000/sin(alf_umb); 
			umbvert = tan(alf_umb)*(y - horiz);

			if vert < umbvert
				shade = 2;
			end
		end
	end
end 

function [eclipse] = eclipseFinding(r, s) 
	theta = mod(asin(6378000/norm(r)), 2*pi);
	psi = mod(acos(dot(r, s)/(norm(r)*norm(s))), 2*pi);

	if psi < theta 
		eclipse = true;
	else 
		eclipse = false;
	end
end