% Function to plot a ground track from latitutes and longitudes - using
% MATLAB's geoplot() and geoscatter() functions. Must be in a new figure or
% on current goegraphic axes.
%
% Author: Griffin Jourda 9/28/22
% 
%	Inputs 
%		lats	:	vector of latitudes (deg) 
%		lons	:	vector of longitudes (deg) 
%
%	Outputs 
%		None

function [] = plotGroundTrack(lats, lons, style) 
	% Plot the orbit ground track - split coordinates to prevent stringing
	for n = 2:length(lons)
		if (sign(lons(n)) == -sign(lons(n - 1))) && abs(lons(n)) > 10
			lons(n) = nan;
		end
	end

	geoplot(lats, lons, style, 'LineWidth', 1)
	hold on 

	% Plot start and end points of the track
	geoscatter(lats(1), lons(1), 'rsquare', 'filled')
	geoscatter(lats(end), lons(end), 'bo', 'filled')












% 	if max(lons) < 175 && min(lons) > -175
% 		geoplot(lats, lons, style, 'LineWidth', 1)
% 	else 
% 		pn = 1;
% 
% 		for n = 2:length(lons)
% 			if sign(lons(n)) == -sign(lons(n-1))
% 				geoplot(lats(pn:n-1), lons(pn:n-1), style, 'LineWidth', 1)
% 				pn = n;
% 			end
% 		end
% 
% 		geoplot(lats(pn:end), lons(pn:end), style, 'LineWidth', 1)
% 	end
end