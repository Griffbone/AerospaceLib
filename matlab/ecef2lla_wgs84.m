% Function to calculate latitude, longitude, and alittude from an ECEF
% position vector using the WGS-84 model. Follows algorithm outlined in:
% https://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf
%
% Author: Griffin Jourda 9/27/2024
%
% Inputs
%	r	:	ECEF position vector (m) 
%
% Outputs
%	lat	:	latitude (rad) 
%	lon	:	longitude (rad)
%	alt	:	ellipsoidal altitude (m)

function [lat, lon, alt] = ecef2lla_wgs84(r)
	% Constants for Earth (WGS-84 model)
	% f = 1/298.257223563;
	% ee = (2*f - f*f);
	a = 6378137.0;
	ee = 6.69437999014e-3;

	% Calculate longitude 
	lon = atan2(r(2), r(1)); 

	% Calcualte prerequisites for lattitude;
	p = sqrt(r(1)^2 + r(2)^2);
	phi_k = atan2(p, r(3));
	delta = 100; 
	
	% Perform loop to find lattitude
	while delta > 1e-9 
		rn = a/sqrt(1 - ee*sin(phi_k)^2);
		h = p/cos(phi_k) - rn;
		phi_kp1 = atan((r(3)/p)*(1 - ee*(rn/(rn + h)))^-1);
		
		delta = abs(phi_kp1 - phi_k);
		phi_k = phi_kp1;
	end
	
	% Compute new height
	alt = p/cos(phi_k) - rn;
	lat = phi_k;
end