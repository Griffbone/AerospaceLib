%% Problem 1
clc; clear; 
load('constants.mat', 'mu')

i = 56;				% inclination, deg
h = 23222e3;		% orbital altitude, m
t = 24;				% number of satellites
p = 3;				% number of planes 
f = 1;				% phasing

[elements, states] = constellation(i, h, t, p, f);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
tspan = linspace(0, 60*60*24, 1000); 

% Simulating Each Satellite 
rs = []; 
for i = 1:t
	r = states(i, 1:3); 
	v = states(i, 4:6); 

	[~, y] = ode45(@(t, y) twoBodyEOM(t, y, mu), tspan, [r'; v'], options);
	rs = [rs, y(:, 1:3)];
end

hold on
for sat = 1:t
	r = rs(:, 1 + 3*(sat-1):3 + 3*(sat-1)); 
	plot3(r(:, 1), r(:, 2), r(:, 3), 'k'); 
	scatter3(r(1, 1), r(1, 2), r(1, 3), 'r', 'filled'); 
end


drawEarth()
xlabel('ECI X [km]')
ylabel('ECI Y [km]')
zlabel('ECI Z [km]')
axis equal 
grid on
view(65, 15)

% Determine Ground contacts
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
	
		if el > 10
			contacts(it, sat) = 1; 
		else 
			contacts(it, sat) = 0; 
		end 
	end
end


figure
tspan = tspan(1:1:end); 
rs = rs(1:1:end, :); 
contacts = contacts(1:1:end, :); 

plot(tspan/60/60, sum(contacts, 2))
xlabel('Time [hrs]')
ylabel('Number of Visible Satellites')
grid on

% M = animateOrbits(tspan, rs, 0.0000001, r_site, contacts);

%% Problem 2 
clc; clear; 
load('constants.mat', 'mu')

i = 97;				% inclination, deg
h = 420;			% orbital altitude, m
t = 100;			% number of satellites
p = 1;				% number of planes 
f = 1;				% phasing

[elements, states] = constellation(i, h, t, p, f);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 

tmax = 1; % tmax (days)
tspan = linspace(0, 60*60*24*tmax, tmax*40e3); 
w = 2*pi/86164;

% Simulating Each Satellite 
lats = [];
lons = [];
rs = []; 
for i = 1:1 % t
	r = states(i, 1:3); 
	v = states(i, 4:6); 

	[~, y] = ode45(@(t, y) twoBodyEOM(t, y, mu), tspan, [r'; v'], options);
	rs = [rs, y(:, 1:3)];

	lati = [];
	loni = [];

	for k = 1:length(tspan)
		r = y(k, 1:3); 
		r_ecef = R3(w*tspan(k))*r'; 
		[lat, lon] = ECEF2LLA(r_ecef', 'gc'); 
		lati = [lati; lat]; 
		loni = [loni; lon];
	end

	lats = [lats, lati]; 
	lons = [lons, loni];
end

% Discretize the earth into 1x1 kilometer sections - only use lattitudes
% +/- 80 degrees. Assume all lines of longitude/lattitude are evenly spaced
% across the whole surface (they all form squares)
h = deg2rad(80*2)*6371;			% height of map +/- 80 deg lat (km)
w = 2*pi*6371;					% width of map +/- 180 deg lon (km)
res = 1;						% grid resolution (km) 
mapgrid = zeros([floor(h/res), floor(w/res)]);
w_cells = size(mapgrid, 2);
h_cells = size(mapgrid, 1);

% Create sub-mask to draw a 24 kilometer diameter circle 
c = drawcircle(24); 

% Cycle through every lat/lon of every satellite and draw a 24 km diameter
% circle around the subsatellite point. Basically, highlight all the cells
% that are imaged by the satellite over the simulation period. 
for sat = 1:1:size(lats, 2)
	lati = lats(:, sat); 
	loni = lons(:, sat); 

	mask = lati <= 80 & lati > -80;
	lati = lati(mask);
	loni = loni(mask); 

	for k = 1:size(lati, 1)
		[i, j] = latlon2cell(lati(k), loni(k), w_cells, h_cells); 

		if ~isnan(i) && ~isnan(j)
			irange = i-12:i+11; 
			jrange = j-12:j+11; 

			if (all(irange) > 0 && all(irange < size(mapgrid, 1))) && ...
					(all(jrange) > 0 && all(jrange < size(mapgrid, 2)))
				mapgrid(i-12:i+11, j-12:j+11) = c + mapgrid(i-12:i+11, j-12:j+11); 
			end

			mapgrid(i, j) = 255; 
		end
	end
end

%% Functions
% Function to return a matrix with a circle with a given diameter inside of
% it (think building a circle in Minecraft)
% 
% Author: Griffin Jourda 12/07/22
% 
%	Inputs 
%		d			:	diameter of the circle 
% 
%	Outptuts 
%		submask		:	dxd matrix of zeros with values of 255 inside of
%						the circle radius 
function [submask] = drawcircle(d)
	s = d; 
	submask = zeros(s); 

	for i = 1:s
		for j = 1:s
			if (i - s/2)^2 + (j - s/2)^2 <= (d/2)^2
				submask(i, j) = 255; 
			end
		end 
	end 
end

% Function to convert lattitude and longitude into indices of an array
% representing discretized portions of the Earth's surface. Hard coded to
% be the subsection between +/- 80 degrees of lattitude. 
% 
% Author: Griffin Jourda 12/7/22
% 
%	Inputs 
%		lat		:	lattitude (deg) 
%		lon		:	longitude (deg) 
%		w		:	width of the earth surface array
%		h		:	height of the earth surface array 
% 
%	Outputs 
%		i		:	vertical index on earth surface array 
%		j		:	horizontal index on earth surface array 

function [i, j] = latlon2cell(lat, lon, w, h) 
	cplon = w/360;			% cells per degree of longitude 
	cplat = h/(80*2);		% cells per degree of lattitude

	j = floor(w/2 + lon*cplon);
	i = floor(h/2 - lat*cplat);

	if i <= 0 || i > h
		i = nan; 
	end

	if j <= 0 || j > w
		j = nan; 
	end
end

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

% Function to plot Earth and axes of an inertial frame on a 3D axes in
% kilometers 
% 
% Author: Griffin Jourda 12/07/22
function [] = drawEarth() 
	[xx, yy, zz] = sphere(20); 
	xx = xx*6371e3; 
	yy = yy*6371e3; 
	zz = zz*6371e3; 

	hold on
	plot3([0, 6371e3*1.5], [0, 0], [0, 0], 'r', 'LineWidth', 2)
	plot3([0, 0], [0, 6371e3*1.5], [0, 0], 'g', 'LineWidth', 2)
	plot3([0, 0], [0, 0], [0, 6371e3*1.5], 'b', 'LineWidth', 2)
	surf(xx, yy, zz, 'FaceAlpha', 0)
	hold off
end