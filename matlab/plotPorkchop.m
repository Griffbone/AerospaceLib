% Function to plot a porkchop plot nicely from previously generated data 
% 
% Author: Griffin Jourda 10/13/22
% 
%	Inputs 
%		xx		:	meshgrid of departure Julian dates 
%		yy		:	meshgrid of arrival Julian dates
%		zz		:	meshgrid of porkchop quantity (C3, excess vel, etc.)
%		levels	:	levels of the quantity to plot 

function [] = plotPorkchop(xx, yy, zz, levels)
	contour(xx, yy, zz, levels)
	colorbar('WestOutside')
	grid on
	hold on
	axis square
	xlabel('Departure Date')
	ylabel('Arrival Date')
	
	% Departure and arrival Julian dates
	x = xx(end,:); 
	y = yy(:, 1); 
	xlim([min(x), max(x)]); 
	ylim([min(y), max(y)]);
	
	% Plot times of flight
	tofmax = max(y) - max(x);
	tofmin = min(y) - max(x);
	tofs = 50:50:tofmax;
	
	for tof = tofs
		plot(x, x + tof, 'k')
		text(x(end), x(end) + tof*1.005, sprintf(' %.3i', tof))
	end
	
	% Conver JDs to calendar dates and do a bunch of tick formatting
	x_ticks = x([1, floor(length(x)/9):floor(length(x)/9):length(x) - floor(length(x)/9), length(x)]);
	y_ticks = y([floor(length(y)/9):floor(length(y)/9):length(y) - floor(length(y)/9), length(y)]);
	
	x_tick_labels = []; 
	y_tick_labels = [];
	
	for i = 1:length(x_ticks)
		[y, m, d, ~, ~, ~, ~] = jd2greg(x_ticks(i));
		str = sprintf('%.2i/%.2i/%.2i', m, d, y - 2000);
		x_tick_labels = [x_tick_labels, {str}];
	end
	
	for i = 1:length(y_ticks)
		[y, m, d, ~, ~, ~, ~] = jd2greg(y_ticks(i));
		str = sprintf('%.2i/%.2i/%.2i', m, d, y - 2000);
		y_tick_labels = [y_tick_labels, {str}];
	end
	
	xticks(x_ticks) 
	yticks(y_ticks)
	yticklabels(y_tick_labels)
	xticklabels(x_tick_labels)
end