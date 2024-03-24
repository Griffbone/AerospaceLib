% Function to retreive orbital elements from a two-line element 
% 
% Author: Griffin Jourda 9/24/22
% 
%	Inputs
%		l1	:	first line of the TLE (string) 
%		l2	:	second line of the TLE (string) 
%	Outputs
%		jd		:	Julian date (day number + day fraction)
%		a		:	semimajor axis (m)
%		e		:	eccentricity 
%		i		:	inclination (deg) 
%		raan	:	right ascension of the ascending node (deg) 
%		w		:	argument of periapsis (deg) 
%		ta		:	true anomaly (deg)
%		epoch	:	string reprsenting epoch and time
%		nd		:	mean motion derivative (rev/day)
%		ndd		:	mean motion derivative (rev/day^2)

function [jd, a, e, i, raan, w, ta, epoch, id, nd, ndd] = readTLE(l1, l2)
	l1 = strsplit(l1, ' ');
	l2 = strsplit(l2, ' ');
	
	id = l1{2};

	% Determine epoch
	yr = eval(l1{4}(1:2));
	[y, ~, ~] = ymd(datetime('today'));
	
	if yr > (y - floor(y/100)*100)
		yr = 1900 + yr; 
	else
		yr = 2000 + yr;
	end
	
	d = eval(l1{4}(3:end));
	jd = greg2jd(yr, 1, 1, 0, 0, 0) + d - 1;

	% Inclination 
	i = eval(l2{3});

	% Right ascension of the ascending node 
	raan = eval(l2{4});

	% Eccentricity 
	e = eval(l2{5})/1e7; 

	% Argument of periapsis
	w = eval(l2{6});

	% Mean anomaly 
	M = eval(l2{7});
	[~, ta] = keplerE(deg2rad(M), e); 
	ta = rad2deg(ta); 

	% Mean motion 
	n = eval(l2{8});
	n = n*(2*pi)*(1/86400);
	
	nd = eval(l1{5}); 
	ndd = eval(l1{6});

	% semimajor axis 
	load('constants.mat', 'mu')
	a = (mu/n^2)^(1/3); % mu from https://ssd.jpl.nasa.gov/astro_par.html 
	a = a;
	
	% epoch string 
	[~, ~, ~, ~, ~, ~, epoch] = jd2greg(jd);
end