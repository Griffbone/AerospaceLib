% Function to get approximate state of a planet from Julian date. Follows
% formulations in: vadimchazov.narod.ru/text_pdf/XSChap8.pdf Table 8.10.2
% and: ssd.jpl.nasa.gov/planets/approx_pos.html.
% 
% Author: Griffin Jourda 10/12/22
% 
%	Inputs 
%		jd	:	Julian date 
%		pn	:	Planet number (1 = Mercury, 2 = Venus, ... 8 = Neptune)
%	
%	Outputs
%		r	:	Position in heliocentric J2000 frame (km) 
%		v	:	Velocity in heliocentric J2000 frame (km/s)
%		mu	:	Planet gravitational parameter (km^3/s^2)

function [r, v, mu] = planetState(jd, pn)
	T = (jd - 2451545)/36525;
	
	% Elements and secular rates for all planets
	switch pn 
		case 1 
			elements = [0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593];
			rates = [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081];
			mu = 22031.868551;
		case 2 
			elements = [0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255]; 
			rates = [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418];
			mu = 324858.592000;
		case 3 
			elements = [1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0];
			rates = [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0];
			mu = 398600.435507;
		case 4 
			elements = [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891];
			rates = [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343];
			mu = 42828.375816;
		case 5 
			elements = [5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909];
			rates = [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106];
			mu = 126712764.100000;
		case 6 
			elements = [9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448];
			rates = [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794];
			mu = 37940584.841800;
		case 7 
			elements = [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503];
			rates = [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589];
			mu = 5794556.400000;
		case 8 
			elements = [30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574]; 
			rates = [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664];
			mu = 6836527.100580;
	end
	
	% Elements at epoch from secular rates
	a = (elements(1) + rates(1)*T)*149597870.691;
	e = (elements(2) + rates(2)*T);
	I = mod(elements(3) + rates(3)*T, 360); 
	L = mod(elements(4) + rates(4)*T, 360); 
	wbar = mod(elements(5) + rates(5)*T, 360); 
	Omega = mod(elements(6) + rates(6)*T, 360); 
	w = wbar - Omega; 
	
	% Mean and true anomaly at epoch 
	M = deg2rad(mod(L - wbar, 360));
	[~, ta] = keplerE(M, e);
	
	% Position and velocity at epoch
	[r, v] = COE2RV(a, e, deg2rad(I), deg2rad(Omega), deg2rad(w), ta, 1.32712440018e11);
end