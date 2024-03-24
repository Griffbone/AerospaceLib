clc; clear; 
load('moontimes.mat', 't_rises', 't_sets')
EOPData = loadEOPData('EOPdata.txt');
lat_station = 33.7488; 
lon_station = -84.3877;

jd0 = min([min(t_rises), min(t_sets)]);
times = zeros(1, length(t_rises) + length(t_sets));

leglabels = [];

if ismember(jd0, t_rises)
	times(1:2:end) = t_rises;
	times(2:2:end) = t_sets;
else
	times(1:2:end) = t_sets;
	times(2:2:end) = t_rises;
end

for i = 1:2:length(times)-1
	t0 = times(i); 
	tf = times(i+1);
	tspan = linspace(t0, tf, 1000);
	
	az_span = []; 
	el_span = [];

	for t = tspan
		[az, el] = moonAzEl(t, lat_station, lon_station, EOPData);
		az_span = [az_span, az];
		el_span = [el_span, el];
	end

	[~, ~, ~, ~, ~, ~, risestr] = jd2greg(t0 - 4/25);
	[~, ~, ~, ~, ~, ~, setstr] = jd2greg(tf - 4/24);

	risestr = formatLegendDate(risestr);
	setstr = formatLegendDate(setstr);

	leglabels = [leglabels, {sprintf('%s | %s', risestr, setstr)}];

	polarplot(az_span*pi/180, el_span, 'LineWidth', 2);
	hold on
end

ax = gca; 
ax.ThetaDir = 'clockwise'; 
ax.ThetaZeroLocation = 'top';
ax.RDir = 'reverse'; 
rlim([0, 90])

rticks([0, 30, 60, 90])
thetaticks([0, 45, 90, 135, 180, 225, 270, 315])
thetaticklabels({'N', '45', 'E', '135', 'S', '225', 'W', '315'});
legend(leglabels)

function [fmt] = formatLegendDate(datestr)
	fmt = [datestr(1:6), datestr(12:17)];
end