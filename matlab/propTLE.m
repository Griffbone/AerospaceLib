% Function to propagate a TLE between given Julian dates using two-body
% motion. Kepler propagator (KPP) option is a little messed up - needs some
% work to match TLEs better over longer periods.
%
% Author: Griffin Jourda 9/28/22
% 
%	Inputs 
%		l1		:	first line of the TLE
%		l2		:	second line of the TLE
%		jd1		:	Julian date of start of propagation period
%		jd2		:	Julian date of end of propagation period 
%		method	:	propagation method ('2B', 'J2', 'KPP')
%
%	Outputs
%		t		:	vector of queried time (s)
%		y_eci	:	array of state vectors in ECI frame (m, m/s)
%		y_ecef	:	array of position vector in ECEF frame (m)
%		lat		:	vector of geocentric satellite latitudes (deg) 
%		lon		:	vector of geocentric satellite longitudes (deg)

function [t, y_eci, y_ecef, lat, lon] = propTLE(l1, l2, jd1, jd2, method)
	load('constants.mat', 'mu')
	options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60);
	
	% Find state vector at epoch
	[jd_epoch, a, e, i, raan, w, ta, ~, ~, nd, ndd] = readTLE(l1, l2);
	[r_epoch, v_epoch] = COE2RV(a, e, deg2rad(i), deg2rad(raan), deg2rad(w), deg2rad(ta), mu);
	
	% Propagate to start of propagation period 
	if jd1 == jd_epoch
		r_0 = r_epoch; 
		v_0 = v_epoch; 
	else		
		% Create propagation time span 
		t = (jd1 - jd_epoch)*86400;
		
		% Propagate to window start 
		if strcmpi(method, '2B') 
			[~, y_epoch2start] = ode45('twoBodyEOM', t, [r_epoch'; v_epoch'], options);
		elseif strcmpi(method, 'J2')
			[~, y_epoch2start] = ode45('twoBodyEOM_J2', t, [r_epoch'; v_epoch'], options);
		elseif strcmpi(method, 'KPP') 
			[~, y_epoch2start, ~] = keplerPropPerturbed(r_epoch, v_epoch, t, 0, 0, mu); 
		end
		
		r_0 = y_epoch2start(end, 1:3); 
		v_0 = y_epoch2start(end, 4:6);
	end
	
	% Propagate to end of propagation period
	t = linspace(0, (jd2 - jd1)*86400, (86400/(90*60)*1000));
	
	if strcmpi(method, '2B') 
		[~, y_out] = ode45('twoBodyEOM', t, [r_0'; v_0'], options);
	elseif strcmpi(method, 'J2')
		[~, y_out] = ode45('twoBodyEOM_J2', t, [r_0'; v_0'], options);
	elseif strcmpi(method, 'KPP') 
		nd = (2*nd)*(2*pi)*(1/86400^2);
		ndd = (6*ndd)*(2*pi)*(1/86400^3);
		[~, y_out, ~] = keplerPropPerturbed(r_0, v_0, t, nd, ndd, mu); 
	end
	
	% Determine ECEF positions and geocentric longitude/latitudes
	y_ecef = []; 
	lats = [];
	lons = [];
	for n = 1:size(y_out, 1) 
		gha = jd2gmst(jd1 + t(n)/86400);
		r_ecef = (R3(gha)*y_out(n, 1:3)')';
		[lat, lon, ~] = ecef2lla_ITRS(r_ecef);
	
		y_ecef = [y_ecef; r_ecef];
		lats = [lats, lat];
		lons = [lons, lon];
	end

	lat = lats; 
	lon = lons;
	y_eci = y_out;
end
