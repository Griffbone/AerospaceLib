%% Functions 
% Function to convert two lat/lon bounds of a swath into a mesh 
function [vertices, tris] = swathTris(lat1, lon1, lat2, lon2)
	vertices = [lat1', lon1'; lat2', lon2'];

	p1 = 1:1:length(lat1) - 1; 
	p2 = 2:1:length(lat1);
	p3 = (1:1:length(lat1) - 1) + length(lat1); 

	p22 = (2:length(lat1)) + length(lat1); 

	tris = [p1', p2', p3'];
	tris = [tris; p3', p22', p2'];
end

% Function to calculate terminal coordinates from initial geographic
% position, heading, and distance 
% 
% Author: Griffin Jourda 12/07/22 
% 
%	Inputs 
%		lat1	:	initial latitude (deg)
%		lon1	:	initial longitude (deg)
%		theta	:	heating (deg)
%		dist	:	distance (km) 
% 
%	Outputs 
%		lat2	:	final latitude (deg)
%		lon2	:	final longtiude (deg)
function [lat2, lon2] = terminalCoords(lat1, lon1, theta, dist)
	lat1 = deg2rad(lat1); 
	theta = deg2rad(theta); 
	del = dist./6371;

	lat2 = asind(sin(lat1).*cos(del) + cos(lat1).*sin(del).*cos(theta));
	lon2 = lon1 + atan2d(sin(theta).*sin(del).*cos(lat1), cos(del) - sin(lat1).*sind(lat2)); 
end

% Function to calculate the initial bearing along the great-circle distance
% between two geographic coordinates 
% 
% Author: Griffin Jourda 12/07/22
% 
%	Inputs 
%		lat1	:	initial latitude (deg) 
%		lon1	:	initial longitude (deg)
%		lat2	:	final latitude (deg)
%		lon2	:	final longitude (deg)
% 
%	Outputs 
%		theta	:	initial bearing between sites
function [theta] = bearing(lat1, lon1, lat2, lon2) 
	lat1 = deg2rad(lat1); 
	lat2 = deg2rad(lat2); 
	lon1 = deg2rad(lon1); 
	lon2 = deg2rad(lon2); 

	dlon = lon2 - lon1; 

	theta = atan2d(sin(dlon).*cos(lat2), cos(lat1).*sin(lat2) - ...
		sin(lat1).*cos(lat2).*cos(dlon));
end

