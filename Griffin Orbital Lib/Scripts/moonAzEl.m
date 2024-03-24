function [az, el] = moonAzEl(jd_utc, lat, lon, EOPData) 
	% Fucntion to calculate elevation of the Moon from a given ground station 
	% 
	% Inputs 
	%	jd_utc	:	utc Julian date 
	%	lat		:	latitude of ground station (deg) 
	%	lon		:	longitude of ground station (deg) 
	%	EOPData	:	Earth orientation parameter data
	% Outputs 
	%	el		:	elevaton of the moon (deg) 

	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData); 
	R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);

	s = moonVector(jd_utc);
	s_ecef = R_icrs_itrf*s;
	s_enu = ECEF2ENU_simple(s_ecef, lat, lon);

	[az, el] = ENU2azel(s_enu);
end