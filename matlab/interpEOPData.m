% Function to interpolate Earth orientation parameter data from arrays.
% Data can be loaded in using loadEOPData function on C04 EOP files. 
% 
% Author: Griffin Jourda 3/17/2023
% 
% Inputs
%	jd_utc	:	utc Julian date
%	EOPdata :	array of EOP data [jd, xp, yp, UT1 - UTC, UT1 - TAI] time
%				values in seconds
% 
% Outputs
%	jd_u1	:	UT1 Julian date
%	jd_tai	:	TAI Julian date
%	jd_tt	:	TT Julian date
%	xp		:	x pole location (arcsec)
%	yp		:	y pole location (arcsec)

function [jd_ut1, jd_tai, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPdata)
	jds = EOPdata(:, 1);
	xps = EOPdata(:, 2);
	yps = EOPdata(:, 3);
	dut1 = EOPdata(:, 4); 
	dtai = EOPdata(:, 5);
	
	if jd_utc > max(jds)
		jd_ut1 = jd_utc + dut1(end)/86400; 
		jd_tai = jd_ut1 - dtai(end)/86400; 
		
		xp = xps(end); 
		yp = yps(end); 
	else
		xp = interp1(jds, xps, jd_utc);
		yp = interp1(jds, yps, jd_utc);

		jd_ut1 = interp1(jds, dut1, jd_utc)/86400 + jd_utc;
		jd_tai = jd_ut1 - interp1(jds, dtai, jd_utc)/86400;
	end
	
	jd_tt = jd_tai + 32.184/86400;
end