function [] = animateAttitude(t, y, kpause, plottrace, I) 
	if plottrace == true
		xis = [];
		for n = 1:length(t)
			R = q2R(y(n, 1:end)');
			xis = [xis; (R*[2; 0; 0])'];
		end
	end

	figure
	hold on
	xlim([-5, 5])
	ylim([-5, 5])
	zlim([-5, 5])
	grid on
	axis equal
	view(-55, 16)
	
	dt = mean(diff(t));
	
	for i = 1:1:size(y, 1)
		cla
		
		q = y(i, 1:end)'; 
		R = q2R(q); 
		
		plotSat(q);
		
		plot3([0, 3], [0, 0], [0, 0], 'r', 'LineWidth', 2)
		plot3([0, 0], [0, 3], [0, 0], 'g', 'LineWidth', 2)
		plot3([0, 0], [0, 0], [0, 3], 'b', 'LineWidth', 2)

		hrs = floor(t(i)/3600); 
		mins = floor(t(i)/60 - hrs*60); 
		secs = t(i) - hrs*3600 - mins*60;

		title(sprintf('T = %i:%i:%.3f s, H = %.4f kg*m^2/s', hrs, mins, secs, norm(I*y(i, 1:3)')))
% 		title(sprintf('T = %i:%i:%.3f', hrs, mins, secs))

		if plottrace == true
			plot3(xis(1:i, 1), xis(1:i, 2), xis(1:i, 3))
		end
		
		pause(dt/kpause)
	end
end

% Function to plot a satellite with a given attitude 
% 
% Author: Griffin Jourda 11/23/22
% 
%	Inputs
%		q	:	attitude quaternion [q1; q2; q3; qs]
function [] = plotSat(q)
	% Cube vertices
	vs = [-1, -1, 1;
		-1, -1, -1; 
		1, -1, -1; %3
		1, -1, 1;
		-1, 1, 1; %5
		-1, 1, -1; 
		1, 1, -1; 
		1, 1, 1]; % 8
	
	% cube surfaces
	S = [3, 7, 8, 4;	% + X
		 6, 2, 1, 5;	% - X
		 7, 6, 5, 8;	% + Y
		 1, 2, 3, 4;	% - Y
		 5, 8, 4, 1;	% + Z
		 2, 3, 7, 6];	% - Z
	 
	% Axes
	as = [2, 0, 0;
		  0, 2, 0;
		  0, 0, 2];
	  
	% Rotate vertices and axes
	R = q2R(q); 
	vs = vs*R';
	as = (as*R');
	
	% Plot faces
	for i = 1:size(S, 1) 
		Si = S(i, :); 
		
		fill3(vs(Si, 1), vs(Si, 2), vs(Si, 3), 'k', 'FaceAlpha', 0)
	end
	
	% Plot axes 
	plot3([0, as(1, 1)], [0, as(1, 2)], [0, as(1, 3)], 'r')
	plot3([0, as(2, 1)], [0, as(2, 2)], [0, as(2, 3)], 'g')
	plot3([0, as(3, 1)], [0, as(3, 2)], [0, as(3, 3)], 'b')
end