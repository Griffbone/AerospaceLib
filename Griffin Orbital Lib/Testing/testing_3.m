clc; clear; close all;
load('constants.mat', 'mu_s', 'au');

% %% 
% jdd1 = 2460260; % greg2jd(2023, 1, 1, 0, 0, 0);
% 
% [vd, va, xx, yy] = porkchop(3, 2, jdd1, jdd1 + 10, jdd1 + 100, jdd1 + 4*365.25, 500);
% %% 
% plotPorkchop(xx, yy, vd, linspace(1, 10, 20))

%% Earth-Venus Transfer
jdd = 2460080; 
jda = 2460260;

[r1, v1] = planetState(jdd, 3); 
[r2, v2, mu_v] = planetState(jda, 2); 
[vhd, vha] = lambertUV(r1, r2, (jda - jdd)*86400, mu_s, 1);

%% Venus-Earth Transfer
clc; 

jda2 = jda + 330;
[r3, ~] = planetState(jda2, 3);
[vhd2, ~] = lambertUV(r2, r3, (jda2 - jda)*86400, mu_s, 1); 


t1 = linspace(0, (jda - jdd)*86400, 1000);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60);
[t1, y1] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), t1, [r1'; vhd'], options); 

t2 = linspace(0, (jda2 - jda)*86400, 1000);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60);
[t2, y2] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), t2, [r2'; vhd2'], options); 

t = [t1; t2 + max(t1)];
y = [y1; y2];

plotTransfer(jdd, t, y, [2, 3])

% disp(norm(vhd - v1)) 
% disp(norm(vha - v2))
% disp(norm(vhd2 - v2))

%% Targeting
clc;
vinf1 = vha - v2
vinf2 = vhd2 - v2

% norm(vinf1)
% norm(vinf2)

% [BR, BT, rp, dvp] = BPlane2(vinf1, vinf2, mu_v)

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

