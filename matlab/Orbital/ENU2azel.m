% Function to determine azimuth and elevation from slant range vector in
% the ENU frame. 
% 
% Author: Griffin Jourda 10/04/22
% 
%	Inputs 
%		r_enu	:	slant range vector in ENU frame 
% 
%	Outputs 
%		az		:	azimuth to satellite (deg) 
%		el		:	elevation to satellite (deg) 

function [az, el] = ENU2azel(r_enu) 
	psi = atan2d(r_enu(2), r_enu(1)); 
	az = 90 - psi;
	el = asind(r_enu(3)/norm(r_enu)); 

	az = mod(az, 360); 
end