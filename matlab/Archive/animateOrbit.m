function [] = animateOrbit(t, y, kpause) 
	figure; hold on; grid on; axis equal; 
	view(-55, 16)

	dt = mean(diff(t));
	w = 2*pi/86164;

	for i = 1:size(y, 1)
		cla
		drawEarth(w*t(i))
		plot3(y(1:i, 1), y(1:i, 2), y(1:i, 3))
		r = y(i, 1:3); 
		scatter3(r(1), r(2), r(3), 'filled', 'r')

		ti = t(i); 
		hrs = floor(ti/3600); 
		mins = floor(ti/60 - hrs*60); 
		secs = ti - hrs*3600 - mins*60; 

		title(sprintf('T = %i:%i:%.3f', hrs, mins, secs))
		pause(dt/kpause)
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