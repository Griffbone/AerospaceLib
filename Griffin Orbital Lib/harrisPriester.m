% Function to determine atmospheric density using the Harris-Priester
% atmosphere model
% 
% Author: Griffin Jourda 3/30/2023
% 
% Inputs
%	h		:	satellite ellipsoidal altitude (m) 
%	r_sat	:	geocentric satellite position vector
%	r_sun	:	geocentric sun position vector
%	inc		:	satellite inclination (rad) 
% Outputs 
%	rho		:	local atmospheric density (kg/m^3)

function [rho] = harrisPriester(h, r_sat, r_sun, inc)
	data = [100 497400.0 497400.0;
			120 24900.0 24900.0;
			130 8377.0 8710.0;
			140 3899.0 4059.0;
			150 2122.0 2215.0;
			160 1263.0 1344.0; 
			170 800.8 875.8;
			180 528.3 601.0;
			190 361.7 429.7;
			200 255.7 316.2;
			210 183.9 239.6;
			220 134.1 185.3;
			230 99.49 145.5;
			240 74.88 115.7;
			250 57.09 93.08;
			260 44.03 75.55;
			270 34.30 61.82;
			280 26.97 50.95;
			290 21.39 42.26;
			300 17.08 35.26;
			320 10.99 25.11;
			340 7.214 18.19;
			360 4.824 13.37;
			380 3.274 9.955;
			400 2.249 7.492;
			420 1.558 5.684;
			440 1.091 4.355;
			460 0.7701 3.362;
			480 0.5474 2.612;
			500 0.3916 2.042;
			520 0.2819 1.605;
			540 0.2042 1.267;
			560 0.1488 1.005;
			580 0.1092 0.7997;
			600 0.08070 0.6390;
			620 0.06012 0.5123;
			640 0.04519 0.4121;
			660 0.03430 0.3325;
			680 0.02632 0.2691;
			700 0.02043 0.2185;
			720 0.01607 0.1779;
			740 0.01281 0.1452;
			760 0.01036 0.1190;
			780 0.008496 0.09776;
			800 0.007069 0.08059;
			840 0.004680 0.05741;
			880 0.003200 0.04210;
			920 0.002210 0.03130;
			960 0.001560 0.02360;
			1000 0.001150 0.01810];
	
	h = h/1000;
	
	if h > 1000
		rho = 0;
	else 
		i = find(data(:, 1) <= h, 1, 'last');
		ip1 = i + 1;

		hi = data(i, 1);
		rmi = data(i, 2);
		rMi = data(i, 3);

		hip1 = data(ip1, 1);
		rmip1 = data(ip1, 2);
		rMip1 = data(ip1, 3);

		% Scale heights
		Hm = (hi - hip1)/log(rmip1/rmi);
		HM = (hi - hip1)/log(rMip1/rMi);

		% Interpolation 
		rm = rmi*exp((hi - h)/Hm);
		rM = rMi*exp((hi - h)/HM);

		% Diurnal density variation
		if (inc >= 1.0471975511966) && (inc <= 2.0943951023932)
			n = 6;
		else 
			n = 2;
		end

		eb = R3(-0.523598775598299)*(r_sun/norm(r_sun));
		er = r_sat/norm(r_sat);
		cosnp = (0.5 + 0.5*(dot(er, eb)))^(0.5*n);

		rho = (rm + (rM - rm)*cosnp)*1e-12;
	end
end