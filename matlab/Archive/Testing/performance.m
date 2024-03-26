clc; clear; 

%% ICRS/ITRS Transformations
clc; clear;
n = 50e3; 
jd0 = greg2jd(2023, 1, 1, 0, 0, 0);
EOPData = loadEOPData('EOPdata.txt');

tic()
for t = 1:n
	jd_utc = jd0 + t/86400;
	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);	% Interpolate EOP data
	R_icrf_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);				% Find ICRF-ITRF matrix
	R_itrf_icrf = R_icrf_itrf';										% Find ITRF-ICRF matrix
end
tf = toc();

fprintf('ICRF-ITRF Transformation Call Time: %.6f ms \n', (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% ICRS/ITRS Transformations with frozen P, N, Pi
clc; clear; 
n = 50e3;
jd_utc = greg2jd(2023, 1, 1, 0, 0, 0);
EOPData = loadEOPData('EOPdata.txt');

[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
T_tt = (jd_tt - 2451545)/36525;
T_tt2 = T_tt^2; 
T_tt3 = T_tt^3; 
D2R = pi/180;

% Calculate precession, nutation, polar motion
zeta = (2306.2181/3600)*T_tt + (0.30188/3600)*T_tt2 + (0.017998/3600)*T_tt3;
nu = (2004.3109/3600)*T_tt - (0.42665/3600)*T_tt2 - (0.041833/3600)*T_tt3;
z = zeta + (0.79280/3600)*T_tt2 + (0.000205/3600)*T_tt3;
[eps, dpsi, deps] = nutationIAU1980(jd_tt);
xp = (xp/3600)*D2R;
yp = (yp/3600)*D2R;

N = R1((-eps - deps))*R3((-dpsi))*R1(eps);
P = R3(-z*D2R)*R2(nu*D2R)*R3(-zeta*D2R);
Pi = [1 0 xp; 0 1 -yp; -xp yp 1];

tic()
for t = 1:n
	GMST = jd2gmst(jd_ut1 + t/86400);
	GAST = GMST + dpsi*cos(eps);
	Theta = R3(GAST);

	R_icrf_itrf = Pi*Theta*N*P;
end
tf = toc();

fprintf('ICRF-ITRF Transformation With Fixed Pi/N/P Call Time: %.6f ms \n', (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% ITRS/LLA Transformations
clc; clear; 
n = 50e3;

rs = randn(n, 3);
rs = rs./vecnorm(rs, 2, 2);
rs = rs*6378e3 + rand(n, 1)*1000e3;

tic()
for i = 1:n
	r = rs(i, :)';
	[lat, lon, alt] = ecef2lla_ITRS(r);
end
tf = toc();

fprintf('ITRF-LLA Transformation Call Time: %.6f ms \n', (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% Nonspherical gravity - MATLAB implementation
clc; clear;
n = 10;
deg = 20;

r = randn(n, 3);
r = (r./vecnorm(r, 2, 2)).*(6378e3 + 500e3 + 100e3.*randn(n, 1));

tic()
for i = 1:size(r, 1)
	ri = r(i, :);
	[gx, gy, gz] = gravitysphericalharmonic(ri, deg);
end
tf = toc();

fprintf('MATAB gravitysphericalharmonic (deg %.i) Call Time: %.6f ms \n', deg, (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% Nonspherical gravity - Library Implementation 
clc; clear;
n = 100;
deg = 20;

r = randn(n, 3);
r = (r./vecnorm(r, 2, 2)).*(6378e3 + 500e3 + 100e3.*randn(n, 1));

tic()
for i = 1:size(r, 1)
	ri = r(i, :);
	[gx, gy, gz] = gravitysphericalharmonic(ri, deg);
end
tf = toc();

fprintf('MATAB gravitysphericalharmonic (deg %.i) Call Time: %.6f ms \n', deg, (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% Sun/Moon position 
clc; clear; 
n = 50e3;
jd0 = greg2jd(2023, 1, 1, 0, 0, 0);

tic()
for t = 1:n
	jd_utc = jd0 + t/86400;
	s_sun = sunVector(jd_utc);
	s_moon = moonVector(jd_utc);
end
tf = toc();

fprintf('Sun and Moon Vector Call Time: %.6f ms \n', (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% Atmosphere
clc; clear;
n = 50e3;

alts = 500e3 + randn(n, 1)*50e3;
incs = rand(n, 1)*pi;
r_sat = randn(n, 3);
r_sat = r_sat./vecnorm(r_sat, 2, 2);

r_sun = randn(n, 3);
r_sun = r_sun./vecnorm(r_sun, 2, 2);

% rho = harrisPriester(926e3, r_sat(1, :)', r_sun(1, :)', incs(1));

tic()
for i = 1:n
	rho = harrisPriester(alts(i), r_sat(i, :)', r_sun(i, :)', incs(i));
end
tf = toc();

fprintf('Harris-Priester Model Call Time: %.6f ms \n', (tf/n)*1000)
fprintf('Total Call Time (%i runs): %.6f s \n', n, tf)

%% High-Fidelity Propagator
clc; clear;

[C, S] = loadSHData('EGM2008.txt', 21);
EOPData = loadEOPData('EOPdata.txt');
rp = 0.63781363e7;
mu = 0.3986004415e15;

%%
t = 44e3;
jd0 = 2459945.5;
beta = 1000/(0.2*1);
y = [5376912.75035258, 1547112.75047136, -3825554.34617017, 4431.295827102, 91.6282171961783, 6254.55997692726]';

tic()
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
toc()*1000