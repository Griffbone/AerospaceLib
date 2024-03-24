clc; clear; close all;
load('constants.mat', 'mu_s', 'au');

%% Transfer planning
BR = -2e5;					% B-Plane B*R target (km)
BT = 0;			% B-Plane B*T target (km)

p1 = 3;		% Departure Planet
p2 = 5;		% Arrival Planet

% Set departure, arrival, and TCM dates
jdd = greg2jd(2022, 6, 0, 0, 0, 0);
jda = jdd + 600;
jd_tcm = jdd + 100;

% Solve Lambert's Problem
[r1, ~] = planetState(jdd, p1); 
[r2, ~] = planetState(jda, p2); 
[v1, ~] = lambertUV(r1, r2, (jda - jdd)*86400, mu_s, 1);

%% Propagation and targeting
% Propagate to TCM maneuver
t1 = linspace(0, (jd_tcm - jdd)*86400, 1000); 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60);
[t1, y1] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), t1, [r1'; v1'], options); 

% Perform B-plane targeting
dv = BPlaneTargeting(jd_tcm, jda, y1(end, 1:3), y1(end, 4:6), BT, BR, p2);
% dv = [0.0143091751931623; -0.00759460635063211; 0.117187754813811];

% Propagate to arrival planet SOI 
t2 = linspace(0, (jda - jd_tcm)*86400, 1000);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60, 'Events', @(t, y) SOIEvent(t, y, jd_tcm, p2));
[t2, y2, te, ~, ~ ] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), t2, [y1(end, 1:3)'; y1(end, 4:6)' + dv], options);

% Propagate about arrival planet
[rp, vp, mu] = planetState(t2(end)/86400 + t1(end)/86400 + jdd, p2);
soi = norm(rp)*(mu/mu_s)^(2/5);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60, 'Events', @(t, y) targetR(t, y, soi));
[t3, y3] = ode45(@(t, y) twoBodyEOM(t, y, mu), [0, 1e9], [y2(end, 1:3)' - rp'; y2(end, 4:6)' - vp'], options);

% % Propagate one period after flyby
% jdd2 = (max(t1) + max(t2) + max(t3))/86400 + jdd;
% [rp, vp, mu] = planetState(jdd2, p2);
% r0 = y3(end, 1:3) + rp; 
% v0 = y3(end, 4:6) + vp; 
% eps = (norm(v0)^2)/2 - mu_s/norm(r0); 
% a = -mu_s/(2*eps); 
% T = 2*pi*sqrt((a^3)/mu);

% t4 = linspace(0, T, 5000);
% options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 86400, 'InitialStep', 60, 'Events', @(t, y) PeEvent(t, y, mu_s));
% [t4, y4] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), t4, [r0'; v0'], options);

%% Plotting

% plotTransfer(jdd, [t1; t2 + max(t1)], [y1; y2], [3, 4])
plotFlyBy(y3, 5)

rps = [];
jda = jdd + (max(t2) + max(t1))/86400;
jds = linspace(jda, jda + max(t3)/86400, size(y3, 1));
y32 = zeros(size(y3));

for i = 1:size(y3, 1)
	[rp, vp] = planetState(jds(i), p2);
	y32(i, 1:3) = y3(i, 1:3) + rp; 
	y32(i, 4:6) = y3(i, 4:6) + vp; 
end
% 
% plotTransfer(jdd, [t1; t2 + max(t1); t3 + max(t2) + max(t1)], [y1; y2; y32], [3, 5])

% 
% 
% 
% % clc
% t = [t1; t2 + max(t1); t3 + max(t2) + max(t1); t4 + max(t3) + max(t2) + max(t1)];
% t = [t1; t2 + max(t1); t4 + max(t3) + max(t2) + max(t1)];
% y = [y1; y2; y3]; %; y4];
% % 
% plotTransfer(jdd, t, y, [p1, p2])

%%
function [] = plotFlyBy(y, pn)
	rps = [2440, 6052, 6371, 3390, 69911, 58232, 25362, 24622];
	rp = rps(pn);
	
	figure
	load('constants.mat', 'mu_s');
	
	% Plot Trajectory
	plot3(y(:,1), y(:,2), y(:,3), 'r-')
	
	hold on 
	axis equal 
	grid on
	
	xlabel('J2000 Ecliptic X [km]')
	ylabel('J2000 Ecliptic Y [km]')
	zlabel('J2000 Ecltipic Z [km]')
	
	scatter3(y(end,1), y(end,2), y(end,3), 'ks')
	
	% Plot planet 
	[xx, yy, zz] = sphere();
	xx = xx*5000;
	yy = yy*5000;
	zz = zz*5000;
	surf(xx, yy, zz, 'FaceColor', 'r', 'FaceAlpha', 0.15, 'EdgeAlpha', 0.25);
	
	% Plot B-plane 
	[rtp, ~, mu] = planetState(0, pn); 
	[~, T, R, B] = BPlaneVecs(y(1,1:3), y(1,4:6), mu);
	appx_soi = norm(rtp)*(mu/mu_s)^(2/5);
	
	rp = 3000; 
	plot3([0, B(1)], [0, B(2)], [0, B(3)], 'k', 'HandleVisibility', 'off')
	plot3([0, T(1)]*rp, [0, T(2)]*rp, [0, T(3)]*rp, 'g', 'HandleVisibility', 'off')
	plot3([0, R(1)]*rp, [0, R(2)]*rp, [0, R(3)]*rp, 'r', 'HandleVisibility', 'off')

	c1 = 2*rp*T + 2*rp*R;
	c2 = 2*rp*T - 2*rp*R; 
	c3 = -2*rp*T - 2*rp*R; 
	c4 = -2*rp*T + 2*rp*R; 
	c5 = c1;
	
	cs = [c1; c2; c3; c4; c5];
	
	fill3(cs(:,1), cs(:,2), cs(:,3), 'k', 'FaceAlpha', 0.1)
	scatter3(B(1), B(2), B(3), 'kx')
	
% 	xlim([-1/3*appx_soi, 1/3*appx_soi])
% 	ylim([-1/3*appx_soi, 1/3*appx_soi])
% 	zlim([-1/3*appx_soi, 1/3*appx_soi])
	
	xlim([-3*rp, 3*rp])
	ylim([-3*rp, 3*rp])
	zlim([-3*rp, 3*rp])
	
	% Legend
	legend({'Trajectory', 'Final Position', 'Planet', 'B-Plane', 'Aymptote B-Plane Intersect'})
end











function [tadot, isterminal, direction] = xyCrossing(~, y)
	tadot = y(3);
	
% 	disp(tadot)
	
	isterminal = 1; 
	direction = 0;
end

% Function to signal an event to MATLAB ode solvers when periapsis is
% reached
% 
% Author: Griffin Jourda 10/14/22
% 
%	Inputs 
%		t	:	integration time (s) 
%		y	:	current state 
%
%	Outputs 
%		tadot		:	dot product of position and eccentricity vector
%		isterminal	:	1 (stop integration when d == 0)
%		direction	:	0 (event can be approached on ascending or
%						descending ta) 

function [tadot, isterminal, direction] = PeEvent(~, y, mu)
	r = y(1:3); 
	v = y(4:6); 
	
	e = (1/mu)*((norm(v)^2 - mu/norm(r))*r - dot(r, v)*v);
	tadot = abs(1 - dot(e, r)/(norm(e)*norm(r)));
	
	if tadot < 1e-6
		tadot = 0; 
% 		disp(acosd(tadot))
	end 
	
	isterminal = 1; 
	direction = 0;
end


% Function to signal an event to MATLAB ode solvers when vehicle is a given
% distance from the central body 
% 
% Author: Griffin Jourda 10/14/22
% 
%	Inputs 
%		t	:	integration time (s)
%		y	:	current state
%		r	:	target distance
%
%	Outputs 
%		del			:	distance from target
%		isterminal	:	1 (stop integration when del == 0) 
%		direction	:	0 (approached from ascending or descending del)

function [del, isterminal, direction] = targetR(~, y, r)
	del = r - norm(y(1:3));
	isterminal = 1; 
	direction = 0;
end

