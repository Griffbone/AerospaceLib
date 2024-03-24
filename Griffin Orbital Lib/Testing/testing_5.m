clc; clear; 
[C, S] = loadSHData('EGM2008.txt', 21);
EOPData = loadEOPData('EOPdata.txt');
rp = 0.63781363e7;
mu = 0.3986004415e15;

%% Run simulation
jd0 = greg2jd(2023, 1, 1, 0, 0, 0);
[r0, v0] = generateSSO(400e3, jd0, 6);
tspan = linspace(0, 86400*1, 1e3);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700);

[t, y] = ode45(@(t, y) EOM_perturbed(t, y, jd0, 1000/(0.2*1), EOPData, C, S, mu, rp), tspan, [r0, v0], options);

%% Back out elements 
elements = zeros(length(t), 6);
sv = [];
for i = 1:length(t)
	r = y(i, 1:3);
	v = y(i, 4:6);
	[a, e, ii, raan, w, ta] = RV2COE(r, v, 398600435436000);
	elements(i, :) = [a, e, ii, raan, w, ta];

	svi = sunVector(t(i)/86400 + jd0);
	sv = [sv; 6378.5e3*svi'/norm(svi)];
end

%% Plotting
tspan2 = linspace(0, max(t), 30*1000);
xx = interp1(t, y(:, 1), tspan2, 'cubic');
yy = interp1(t, y(:, 2), tspan2, 'cubic');
zz = interp1(t, y(:, 3), tspan2, 'cubic');

plotInertialAxes()
plot3(xx/1000, yy/1000, zz/1000)
plot3(sv(:, 1)/1000, sv(:, 2)/1000, sv(:, 3)/1000, 'k', 'LineWidth', 1)

%%

jd0 = greg2jd(2023, 1, 1, 0, 0, 0);
r_sun = sunVector(jd0);
r_sat = (-r_sun/norm(r_sun))*(6378e3 + 500e3)

% eclipse()



%% Functions
% beta - ballistic coefficient (m/CdA)
function [ydot] = EOM_perturbed(t, y, jd0, beta, EOPData, C, S, mu, rp)	
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
	
	% Drag
	h = cross(r_eci, v_eci);
	inc = acos(h(3)/norm(h));
	[~, ~, alt] = ecef2lla_ITRS(r_ecef);
	rho = harrisPriester(alt, r_eci, s_sun, inc);
	a_d = -(1/2)*rho*norm(v_eci)*v_eci*(1/beta);
	
	% Total state derivative 
	a = a_earth + a_sun + a_moon + a_d;
	ydot = [v_eci; a];
end