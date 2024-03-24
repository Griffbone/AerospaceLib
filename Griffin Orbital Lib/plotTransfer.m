% Function to plot a transfer in the solar system 
% 
% Author: Griffin Jourda 10/13/22
% 
%	Inputs
%		j0	:	departure Julian date
%		t	:	array of time from integrator (s) 
%		y	:	spacecraft state array (km, km/s)
%		ps	:	vector of planet numbers (1 = Mercury, 2 = Venus, ... 8 =
%				Neptune)
% 
%	Outputs 
%		

function [] = plotTransfer(j0, t, y, ps)
	load('constants.mat', 'mu_s', 'au');
	cs = ['k', 'k', 'b', 'r', 'k', 'k', 'k', 'k'];
	pnames = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
	jds = t./86400 + j0;
	
	% Plot transfer orbit
	figure
	plot3(y(:,1), y(:,2), y(:,3), 'r', 'LineWidth', 0.75);
	hold on 
	axis equal 
	grid on
	zlim([-au, au])
	
	llabels = [{'Trajectory'}];
	
	% Plot planets 
	for pn = ps
		rps = []; 

		for n = 1:length(jds)
			[rp, ~, mu] = planetState(jds(n), pn); 
			rps = [rps; rp];
		end
		
		plot3(rps(:,1), rps(:,2), rps(:,3), cs(pn), 'LineWidth', 0.75)

		% Plot SOI at end of propagation period 
		soi = norm(rps(end, :))*(mu/mu_s)^(2/5);

		[xx, yy, zz] = sphere(); 
		xx = xx*soi + rps(end, 1); 
		yy = yy*soi + rps(end, 2); 
		zz = zz*soi + rps(end, 3); 

		surf(xx, yy, zz, 'Facecolor', cs(pn), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.25)
		
		llabels = [llabels, strcat(pnames(pn), ' Orbit'), strcat(pnames(pn), ' SOI')];
	end

	% Plot Sun
	[xx, yy, zz] = sphere();
	xx = xx*696340;
	yy = yy*696340;
	zz = zz*696340;
	surf(xx, yy, zz, 'FaceColor', 'y', 'EdgeAlpha', 0.25);
	
	xlabel('J2000 Ecliptic X [km]')
	ylabel('J2000 Ecliptic Y [km]')
	zlabel('J2000 Ecliptic Z [km]')
	
	llabels = [llabels, 'Sun'];
	legend(llabels)
end

% function [] = plotline(p, style)
% 	plot3([0, p(1)], [0, p(2)], [0, p(3)], style)
% end