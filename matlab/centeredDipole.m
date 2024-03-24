% Function to calculate local magnetic field given ECEF coordinates using
% the centered dipole magnetic model 
% 
% Author: Griffin Jourda 1/17/2023
%
%	Inputs
%		r_ecef	:	3x1 ECEF position vector (m)
%		year	:	decimal year corresponding to epoch of the position
%					vector 
%
%	Outputs 
%		B_ecef	:	magnetic field vector expressed in ECEF coordinates
%					(nT)

function [B_ecef] = centeredDipole(r_ecef, year)
	% Model parameters 
	[lat, lon, m] = dipoleData(year); 
	lat = deg2rad(lat);
	lon = deg2rad(lon); 

	% Convert to geomagnetic coordinates 
	R_ecef_gm = R2(pi/2 - lat)*R3(lon);
	r_gm = R_ecef_gm*r_ecef;

	lon_gm = atan2(r_gm(2), r_gm(1)); 
	lat_gm = atan2(r_gm(3), sqrt(r_gm(1)^2 + r_gm(2)^2));
	r = norm(r_ecef); 

	% Determine field parameters 
	mu0 = (4*pi)*1e-7;		% N/A^2
	k0 = (mu0*m)/(4*pi);	% T*m^3

	% Calculate field vector
	Br = (-2*k0/r^3)*sin(lat_gm); 
	Blam = (k0/r^3)*cos(lat_gm); 
	B_enu_gm = [0; Blam; Br];

	R_topo_gm = [-sin(lon_gm), cos(lon_gm), 0; 
		-sin(lat_gm)*cos(lon_gm), -sin(lat_gm)*sin(lon_gm), cos(lat_gm); 
		cos(lat_gm)*cos(lon_gm), cos(lat_gm)*sin(lon_gm), sin(lat_gm)]';

	B_gm = R_topo_gm*B_enu_gm;

	% Convert back to ECEF 
	B_ecef = R_ecef_gm'*B_gm;
	B_ecef = B_ecef*1e9;		% convert to nT
	B_ecef(2) = -B_ecef(2);		% flip the second element (IDK why, but it works)
end

function [lat, lon, m] = dipoleData(year)
	years = [1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, 1945, ...
		1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, 2000, ... 
		2005, 2010, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, ...
		2024, 2025];

	lats = [78.7, 78.7, 78.7, 78.6, 78.6, 78.6, 78.6, 78.6, 78.5, 78.5, ...
		78.5, 78.5, 78.6, 78.6, 78.7, 78.8, 78.9, 79.0, 79.2, 79.4, 79.6, ...
		79.8, 80.1, 80.4, 80.4, 80.5, 80.5, 80.6, 80.7, 80.7, 80.7, 80.8, ...
		80.8, 80.9];	% north latitude

	lons = [68.8, 68.7, 68.7, 68.6, 68.4, 68.3, 68.3, 68.4, 68.5, 68.5, ...
		68.8, 69.2, 69.5, 69.9, 70.2, 70.5, 70.8, 70.9, 71.1, 71.4, 71.6, ...
		71.8, 72.2, 72.6, 72.6, 72.6, 72.7, 72.7, 72.7, 72.7, 72.7, 72.7, ...
		72.6, 72.6];	% west longitude 

	ms = [8.32, 8.3, 8.27, 8.24, 8.2, 8.16, 8.13, 8.11, 8.09, 8.08, 8.06, ...
		8.05, 8.03, 8, 7.97, 7.94, 7.91, 7.87, 7.84, 7.81, 7.79, 7.77, 7.75, ...
		7.72, 7.72, 7.72, 7.71, 7.71, 7.71, 7.71, 7.7, 7.7, 7.7, 7.7];
	
	lat = interp1(years, lats, year); 
	lon = interp1(years, lons, year); 
	m = interp1(years, ms, year)*1e22; 
end