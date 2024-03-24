clc; clear; 
load('constants.mat', 'mu')

i = 97;				% inclination, deg
h = 400e3;			% orbital altitude, m
t = 100;			% number of satellites
p = 1;				% number of planes 
f = 1;				% phasing

[elements, states] = constellation(i, h, t, p, f);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
tspan = linspace(0, 60*60*24, 1000); 

%% Simulating Each Satellite 

rs = []; 
for i = 1:t
	r = states(i, 1:3); 
	v = states(i, 4:6); 

	[~, y] = ode45(@(t, y) twoBodyEOM(t, y, mu), tspan, [r'; v'], options);
	rs = [rs, y(:, 1:3)];
end

% rs = rs(1, :); 
% tspan = tspan(1); 

%% Ground contacts
lat = 33.77206497715465; lon = -84.39583438822801; 
r_site = LLA2ECEF_simple(lat, lon, 0); 

contacts = [];
els = [];
w = 2*pi/86164; 

for it = 1:length(tspan) 
	for sat = 1:t
		r_eci = rs(it, 1 + 3*(sat - 1):3 + 3*(sat - 1)); 
		r_ecef = R3(w*tspan(it))*r_eci'; 
		r_enu = ECEF2ENU_simple(r_ecef', lat, lon); 

		[~, el] = ENU2azel(r_enu); 
	
		if el > 20 
			contacts(it, sat) = 1; 
		else 
			contacts(it, sat) = 0; 
		end 
	end
end

% plot(tspan/60/60, sum(contacts, 2))
% plot(tspan/60^2, els)
% 
M = animateOrbits(tspan, rs, 0.0000001, r_site, contacts);

%% Functions
% Functiont to calculate the initial orbit elements for a satellite
% constellation using Walker's delta notation 
% 
% Author: Griffin Jourda 12/01/22
% 
%	Inputs 
%		i	:	inclination (deg) 
%		h	:	orbital altitude (m)
%		t	:	total number of satellites 
%		p	:	number of evenly-spaced orbital planes 
%		f	:	phase difference between the planes 
%
%	Outputs 
%		elements	:	tx6 array of orbital elements for each satellite 
%						[a, e, i, raan, w, ta] (a in m, angles in rad)
%		states		:	tx6 array of orbital state for each satellite 
%						[x, y, z, vx, vy, vz] (r in m, v in m/s)
function [elements, states] = constellation(i, h, t, p, f)
	sats_per_plane = t/p; 
	plane_spacing = deg2rad(360/sats_per_plane); 
	node_spacing = deg2rad(360/p);
	phasing = deg2rad(360/t)*f; 

	a = 6371e3 + h; 
	i = deg2rad(i); 
	
	elements = []; 
	states = [];
	for planes = 1:p
		raan = (planes - 1)*node_spacing;
		ta0 = (planes - 1)*phasing; 
		ta = 0; 

		for sats = 1:sats_per_plane
			ta = ta0 + plane_spacing*(sats - 1); 
			
			[r, v] = COE2RV(a, 0, i, raan, 0, ta, 398600435507000);
			elements = [elements; a, 0, i, raan, 0, ta];
			states = [states; r, v]; 
		end
	end
end