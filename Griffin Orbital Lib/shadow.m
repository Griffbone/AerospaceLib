% Function to determine if a satellite is in Earth's shadow
% 
% Author: Griffin Jourda 4/6/2023
% 
% Inputs 
%	s_sun	:	geocentric sun vector 
%	s_Sat	:	geocentric satellite vector 
% Outputs 
%	eclipse	:	logical value denoting eclipse (1 = eclipse, 0 = no eclipse)

function [eclipse] = shadow(s_sun, s_sat)
% 	s = norm(s_sat);
	r_b = 6378e3;
	s0 = -(s_sat'*s_sun)/norm(s_sun);
% 	l = sqrt(s^2 - s0^2);

	nss = norm(s_sun);
% 	sf1 = (r_sun + r_b)/nss; 
	sf2 = (696340000 - 6378140)/nss;

% 	c1 = s0 + r_b/sf1;
	c2 = s0 - r_b/sf2;

% 	l1 = c1*tan(asin(sf1));
	l2 = c2*tan(asin(sf2));

	if l2 < 0
		eclipse = true;
	else 
		eclipse = false;
	end
end