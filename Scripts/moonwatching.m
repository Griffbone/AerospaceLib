clc; clear;
EOPData = loadEOPData('EOPdata.txt');

%% Inputs - Start and End Dates, Station Location
clc;
lat_station = 33.7488; 
lon_station = -84.3877;

utc = datetime('now', 'TimeZone', 'Z');
jd_now = greg2jd(utc.Year, utc.Month, utc.Day, utc.Hour, utc.Minute, utc.Second);
jd_now = greg2jd(2023, 4, 1, utc.Hour, utc.Minute, utc.Second);
[az, el] = moonAzEl(jd_now, lat_station, lon_station, EOPData);

fprintf('Current Moon Azimuth : %.2f deg\n', az)
fprintf('Current Moon Elevaion: %.2f deg\n', el)

tspan = linspace(jd_now, jd_now + 10, 10000);

%% Searching for moon rises
el = []; 
t_rises = [];
t_sets = [];

for jd_utc = tspan
	[~, eli] = moonAzEl(jd_utc, lat_station, lon_station, EOPData);
	el = [el, eli];
end

for i = 2:length(el) 
	if el(i) > 0 && el(i-1) < 0
		t_rises = [t_rises, (tspan(i) + tspan(i-1))/2];
		[~, ~, ~, ~, ~, ~, timestr] = jd2greg(t_rises(end) - 4/24);
		disp(timestr(1:end-6))
	end

	if el(i) < 0 && el(i-1) > 0
		t_sets = [t_sets, (tspan(i) + tspan(i-1))/2];
		[~, ~, ~, ~, ~, ~, timestr] = jd2greg(t_sets(end) - 4/24);
		disp(timestr(1:end-6))
	end
end

save('moontimes.mat', 't_sets', 't_rises');