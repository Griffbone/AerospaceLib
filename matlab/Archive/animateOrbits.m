% Function to animate the orbit of multiple satellites 
% 
% Author: Griffin Jourda 12/02/22 
% 
%	Inputs 
%		t		:	time (s) 
%		y		:	ix3*n matrix of satellite x/y/z positions in meters - i
%					represents the number of states (same length as t) while n
%					repreesents the number of satellites
%		k		:	multiple of time step to pause between plot draws 
%		site	:	ECI coordinates of a ground station 
%		contacts:	ixn matrix of satellite contacts - 1 if the satellite
%					sees the ground station and 0 if it does not 
function [M] = animateOrbits(t, y, k, site, contacts)
	h = figure; hold on; grid on; axis equal; 
	h.Visible = 'off';
	view(-55, 16)

	dt = mean(diff(t));
	w = 2*pi/86164;
	M = [];

	for i = 1:size(y, 1)
		cla
		drawEarth(w*t(i))

		if nargin == 5 
			R_site = R3(-w*t(i))*site';
			scatter3(R_site(1), R_site(2), R_site(3), 'filled', 'k')
		end

		for idx = 1:3:size(y, 2)
			r = y(i, idx:idx+2); 
			rsi = y(1:i, idx:idx+2);
			plot3(rsi(:, 1), rsi(:, 2), rsi(:, 3), 'k')
			scatter3(r(1), r(2), r(3), 'filled', 'r')

			if nargin == 5
				if contacts(i, (idx+2)/3) == 1
					plot3([R_site(1), r(1)], [R_site(2), r(2)], [R_site(3), r(3)], 'k--');
				end
			end 
		end

		ti = t(i); 
		hrs = floor(ti/3600); 
		mins = floor(ti/60 - hrs*60); 
		secs = ti - hrs*3600 - mins*60; 

		title(sprintf('T = %i:%i:%.3f', hrs, mins, secs))

% 		M = [M; getframe()];
% 		imwrite(gcf, 'galileo.gif', 'gif', 'WriteMode', 'append');
		exportgraphics(gcf,'galileo.gif','Append',true);
		waitbar(i/size(y, 1))	
	end
end

function [] = drawEarth(theta) 
	[xx, yy, zz] = sphere(20); 
	xx = xx*6378e3; 
	yy = yy*6378e3; 
	zz = zz*6378e3; 

	R = R3(-theta); 
	y = R*[0; 6378e3; 0]*1.2; 
	x = R*[6378e3; 0; 0]*1.2; 

	plot3([0, x(1)], [0, x(2)], [0, x(3)], 'r', 'LineWidth', 2)
	plot3([0, y(1)], [0, y(2)], [0, y(3)], 'g', 'LineWidth', 2)
	plot3([0, 0], [0, 0], [0, 1.2*6378e3], 'b', 'LineWidth', 2)

	surf(xx, yy, zz, 'FaceAlpha', 0)
end