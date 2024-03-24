% Function to determine whether a satellite in orbit of Earth is in
% Eclipse. Follows algorithm 34 in Vallado. 
%
% Author: Griffin Jourda 3/7/2023
% 
%	Inputs 
%		r		:	ECI satellite position vector 
%		r_sun	:	ECI sun position vector 
%
%	Outputs 
%		r_sun	:	sun vector in ECI frame (m)

function [eclipse] = earthEclipse(r, r_sun)
	% Determine umbra/penumbra shadow sizes
	r_s = 695700000;
	r_e = 6378000;
	R = norm(r_sun); 
% 	alf_umb = asin((r_s - r_e)/R);
	alf_pen = asin((r_s + r_e)/R);

	eclipse = false;
	xi = acosd(-r_sun, r)/(norm(r_sun)*norm(r));
	h = norm(r)*cos(xi);
	v = norm(r)*sin(xi); 
	x = r_e/sin(alf_pen);
	p_v = tan(alf_pen)*(x + h);

	if v <= p_v
		eclipse = true;
	end
end