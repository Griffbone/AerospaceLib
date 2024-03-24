% Function to plot a spherical planet and axes for an inertial frame 
% 
% Author: Griffin Jourda 12/01/22
% 
%	Inputs
%		r	:	planet radius 
%		f	:	axis length as fraction of planet radius 
function [] = plotInertialFrame(r, f)
	[xx, yy, zz] = sphere(25); 
	xx = xx*r; 
	yy = yy*r; 
	zz = zz*r; 

	hold on;
	surf(xx, yy, zz, 'FaceColor', 'none');
	plot3([0, r*f], [0, 0], [0, 0], 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off')
	plot3([0, 0], [0, r*f], [0, 0], 'g', 'LineWidth', 1.5, 'HandleVisibility', 'off')
	plot3([0, 0], [0, 0], [0, f*r], 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off')
end 