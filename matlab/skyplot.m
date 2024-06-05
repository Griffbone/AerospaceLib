% Function to plot a skyplot given azimuth and elevation 
% 
% Author: Griffin Jourda 10/06/22
% 
%	Inputs 
%		az	:	azimuth (deg) 
%		el	:	elevation (deg) 
% 
%	Outputs 
%		

function [] = skyplot(az, el)
	figure

	% idx = find(diff(az) == max(diff(az)));
	% az(idx) = nan; 
	% el(idx) = nan; 

	polarscatter(deg2rad(az), el, 'r', 'filled')

	ax = gca; 
	ax.ThetaDir = 'clockwise'; 
	ax.ThetaZeroLocation = 'top';
	ax.RDir = 'reverse'; 
	rlim([0, 90])

	rticks([0, 30, 60, 90])
	thetaticks([0, 45, 90, 135, 180, 225, 270, 315])
	thetaticklabels({'N', '45', 'E', '135', 'S', '225', 'W', '315'});
end

