function plotInertialAxes(fname)
	[xx, yy, zz] = ellipsoid(0, 0, 0, 6378130/1000, 6378130/1000, 6371009/1000, 30);
	surf(xx, yy, zz, 'FaceColor', 'w', 'EdgeColor', [0.5 0.5 0.5], 'HandleVisibility', 'off');
	hold on
	axis equal
	
	plot3([0 1]*6378130/1000*1.5, [0 0], [0 0], 'r', 'LineWidth', 2, 'HandleVisibility', 'off')
	plot3([0 0], [0 1]*6378130/1000*1.5, [0 0], 'g', 'LineWidth', 2, 'HandleVisibility', 'off')
	plot3([0 0], [0 0], [0 1]*6378130/1000*1.5, 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
	
	xlabel('X [km]')
	ylabel('Y [km]')
	zlabel('Z [km]')
	
	if nargin ~= 0
		data = readmatrix(fname);
		data = unique(data, 'rows');
		
		border_ecef = zeros(size(data, 1), 3);

		for i = 1:size(data, 1) 
			lon = data(i, 1);
			lat = data(i, 2); 

			r = LLA2ECEF(lat, lon, 0);
			border_ecef(i, :) = r;	
		end
		
		plot3(border_ecef(:, 1)/1000, border_ecef(:, 2)/1000, border_ecef(:, 3)/1000, '.k', 'MarkerSize', 2)
		axis equal
	end 
end