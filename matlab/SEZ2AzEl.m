% Function to convert SEZ coordiantes to azimuth and elevation angles
% 
% Author: Griffin Jourda 10/01/22
% 
%	Inputs
%		r_sez	:	SEZ position vector of the satellite
%
%	Outputs
%		az		:	satellite azimuth (deg) 
%		el		:	satellite elevation (deg) 

function [az, el] = SEZ2AzEl(r_sez) 
	el = atan2d(r_sez(3), sqrt(r_sez(1)^2 + r_sez(2)^2));
	
	sb = r_sez(2)/(sqrt(r_sez(1)^2 + r_sez(2)^2));
	cb = -r_sez(1)/(sqrt(r_sez(1)^2 + r_sez(2)^2));	
	az = atan2d(sb, cb);
	
	az = mod(az, 360);
end