clc; clear; 
load('constants.mat', 'mu'); 
r0 = [6378e3, 0, 0];
v0 = [0, 8000, 500];
y0 = [r0'; v0'];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
[t, y] = ode45(@(t, y) twoBodyEOM(t, y, mu), linspace(0, 10e3, 1e3), y0, options); 

tic()
% M = animateOrbit(t, y);
M = animateOrbit(t, y(:, 1:3));
toc()

%% 
figure
axis off 
movie(M)

function [M] = animateOrbit(t, rs) 
	f = figure;
% 	ax = axes();
	h.Visible = 'off';
% 	ax.Visible = 'off';
	view(-55, 16); 
	
	% Draw earth
	[xx, yy, zz] = sphere(20); 
	xx = xx*6378e3; 
	yy = yy*6378e3; 
	zz = zz*6378e3;
	surf(xx, yy, zz, 'FaceAlpha', 0)

	% Create a point object for each satellite
	ps = []; 
	for i = 1:3:size(rs, 2) 
		r = rs(1, i:i+2); 
		p = scatter3(r(1), r(2), r(3), 'r'); 
		ps = [ps, p]; 
	end

	% Loop through time series and update satellite points 
	M(length(t)-1) = struct('cdata', [], 'colormap', []);

	for i = 2:length(t)
		disp(i/length(t)); 

		xlim([-50e3, 50e3])
		ylim([-50e3, 50e3])
		zlim([-50e3, 50e3])

		for j = 1:length(ps)
			r =  rs(i, 1 + 3*(j-1):3 + 3*(j-1));
			p = ps(j); 
			p.XData = r(1); 
			p.YData = r(2); 
			p.ZData = r(3); 
		end
		
		M(i-1) = getframe; 
	end 
end