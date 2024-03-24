clc; clear; close all; 

%% Propagate tLE
l1 = '1 49260U 21088A   22299.20586267  .00000286  00000-0  73615-4 0  9992';
l2 = '2 49260  98.2150   7.7454 0000959 104.0799 256.0506 14.57114149 57315';
jd = readTLE(l1, l2);

[jd, ~, ~, ~, ~, ~, ~, ~] = readTLE(l1, l2);
[t, y_eci, y_ecef, ~, ~] = propTLE(l1, l2, jd, jd + 1, 'J2');

%% Ground pass analysis
station_coords = [43.5460, -96.7313, 0; 
				  64.9777, -147.5, 0;
				  77.8750, 20.9752, 0;
				  -23.6980, 133.8807, 0;
				  53.3601, 13.0730, 0];

lam = station_coords(4, 2);
longVL = lam; 
phi = station_coords(4, 1);
latVL = phi;

% r_e = 6371000;
% ECEFvl = [r_e*sind(latVL)*cosd(longVL); r_e*sind(latVL)*sind(longVL); r_e*cosd(latVL)];

R = [-sind(lam), cosd(lam), 0; 
	-sind(phi)*cosd(lam), -sind(phi)*sind(lam), cosd(phi); 
	cosd(phi)*cosd(lam), cosd(phi)*sind(lam), sind(phi)];

al = [];
ep = []; 

for i = 1:size(y_ecef, 1) 
% 	r_ecef = y_ecef(i, :);
% 	r_enu = ECEF2ENU_simple(r_ecef, phi, lam);

% 	r_enu = ECEF2ENU_simple(y_ecef(i, :), phi, lam);
	r_enu = R*y_ecef(i, :)' - [0; 0; 6378e3];
	
	[a, e] = ENU2azel(r_enu);
	al(i) = a; 
	ep(i) = e;
end

% for i = 1:1:size(LOStopo, 2)
%     al(i) = 90 - atan2d(LOStopo(2,i), LOStopo(1,i));
%     ep(i) = asind(LOStopo(3,i)/norm(LOStopo(:, i)));
% end % now we have a bunch of LOS vector

mask = ep > 0;

skyplot(al(mask), ep(mask))

%%
% elevMask = ep >= (5);
% mask2 = LOStopo(3,:) > 0;
% mask2 = logical(elevMask .* mask2);
% plotEp = ep(mask2);
% plotAl = al(mask2);
% timesVisible = t(mask2);
% figure
% polarplot(deg2rad(plotAl), plotEp, '.')
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% ax.RDir = 'reverse';
% rlim([0 90]);
% rtickangle(10)
% thetaticks([0 45 90 135 180 225 270 315])
% thetaticklabels({'North', '', 'East', '', 'South', '', 'West', ''})
% title("Orbit Skyplot", GSname)