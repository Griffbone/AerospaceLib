% Wrapper function to run Vallado's SGP4 implementation for a given TLE
% 
% Author: Griffin Jourda 4/07/2023
%
% Inputs 
%	tle			:	filename of properly formatted TLE
%	MJD_Epoch2	:	MJD of epoch to start the propagation of the satellite
%					(propagate to this date, then return data for tspan after it)
%	tspan		:	span of seconds to return data for after being
%					propagated to MJD_Epoch2
% Outputs
%	r_eci		:	ECI position vectors (km)
%	v_eci		:	ECI velocity vectors (km/s)
%	r_ecef		:	ECEF position vectors (km)
%	v_ecef		:	ECEF velocity vectors (km/s)
%	MJD_EPOCH	:	Modified julian date of TLE epoch

function [r_eci, v_eci, r_ecef, v_ecef, MJD_Epoch] = runSGP4(tle, MJD_Epoch2, tspan)
	addpath('SGP4');
	global const % SAT_Const

	% set up satellite and determine epoch
	[satdata, year, doy]  = loadTLE(tle);

	if (year < 57)
		year = year + 2000;
	else
		year = year + 1900;
	end

	[mon,day,hr,minute,sec] = days2mdh(year,doy);
	MJD_Epoch = Mjday(year,mon,day,hr,minute,sec);
		
	if isnan(MJD_Epoch2)
		tstart = 0;
	else 
		tstart = (MJD_Epoch2 - MJD_Epoch)*1440;
	end

	% read eop data 
	fid = fopen('eop19620101.txt','r');
	eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
	fclose(fid);

	% run sgp4
	y_eci = zeros(length(tspan), 6);
	y_ecef = zeros(length(tspan), 6);
	
	for i = 1:length(tspan)
		tsince = tstart + tspan(i)/60;
		[rteme, vteme] = sgp4(tsince, satdata);

		% time conversions
		MJD_UTC = MJD_Epoch+tsince/1440;

		% eop interpolation
		[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,~,~,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
		[~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
		MJD_UT1 = MJD_UTC + UT1_UTC/86400;
		MJD_TT  = MJD_UTC + TT_UTC/86400;
		T = (MJD_TT-const.MJD_J2000)/36525;

		% frame transformations
		[reci, veci] = teme2eci(rteme,vteme,T,dpsi,deps);
		[recef,vecef] = teme2ecef(rteme,vteme,T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2);

		y_eci(i, :) = [reci', veci'];
		y_ecef(i, :) = [recef', vecef'];
	end
	
	r_eci = y_eci(:, 1:3);
	v_eci = y_eci(:, 4:6);
	r_ecef = y_ecef(:, 1:3);
	v_ecef = y_ecef(:, 4:6);
end
