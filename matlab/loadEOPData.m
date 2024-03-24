% Function to load Earth Orientation parameter data from C04 files retrived
% from: https://hpiers.obspm.fr/eop-pc/index.php?index=C04&lang=en
% 
% Author: Griffin Jourda 3/17/2023
% 
% Inputs
%	fname	:	file storing EOP data
% 
% Outputs 
%	data	:	nx5 array of EOP data [jd, xp (arcsec), yp (arcsec), UT1 -
%				UTC (sec), UT1 - TAI (s)]

function [data] = loadEOPData(fname)
	% Read in data from file
	data = readmatrix(fname); 

	% Parse data
	jd = greg2jd(data(:, 1), data(:, 2), data(:, 3), 0, 0, 0);
	xp = data(:, 4)/1000;
	yp = data(:, 6)/1000;
	dut1 = data(:, 8)/1000;		% UT1 - UTC = dut1
	dtai = data(:, 10)/1000;	% UT1 - TAI = dtai

	% Store in array 
	data = [jd, xp, yp, dut1, dtai];
end