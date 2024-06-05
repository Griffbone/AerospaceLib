clc; clear; 
addpath('..\..\matlab\')

%% Propagate each TLE to date 
jd0 = greg2jd(2024, 7, 13, 9, 0, 0);
jdf = greg2jd(2024, 7, 13, 9 + 7, 0, 0);

r_ends = [];
prns = [];
satnums = [];

for i = 1:31
	fid = fopen(['tles\tle_', num2str(i), '.txt']);
	l0 = fgetl(fid);
	l1 = fgetl(fid);
	l2 = fgetl(fid);
	fclose(fid);

	prns = [prns, str2num(regexp(l0, '(?<=PRN\s)\d{2}', 'match', 'once'))];

	[~, y_eci, y_ecef, ~, ~] = propTLE(l1, l2, jd0, jdf, 'J2');
	r_ends = [r_ends; y_ecef(end, 1:3)];
end 

%% 
lat_gs = 35.347566370233956*pi/180;
lon_gs = -117.8094102639529*pi/180;
alt_gs = 1068.0192;

lat_gt = 33.772348601776734*pi/180;
lon_gt = -84.39470938244997*pi/180;
alt_gt = 297;

r_gs = LLA2ECEF(lat_gs, lon_gs, alt_gs);
r_gt = LLA2ECEF(lat_gt, lon_gt, alt_gt);
r_ends = [r_ends; r_gt];
prns = [prns, 100];
plotInertialAxes(); 

azelr = [];
prn_act = [];

for i = 1:length(prns)
	r = r_ends(i, :);
	s = r - r_gs;
	
	el = pi/2 - acos(dot(r_gs/norm(r_gs), s/norm(s)));
    r_enu = ECEF2ENU(r', lat_gs, lon_gs, alt_gs);
    [az, el2] = ENU2azel(r_enu);

    % disp([el - el2]*180/pi)
    
    disp(prns(i))
	if (el2 > 0) || (prns(i) == 100)
		plot3([r_gs(1), r(1)]/1000, [r_gs(2), r(2)]/1000, [r_gs(3), r(3)]/1000, 'k')
		scatter3(r(1)/1000, r(2)/1000, r(3)/1000, 'k', 'filled', 'square');
		% text(r(1)/1000, r(2)/1000, r(3)/1000, num2str(prns(i)), 'color', 'red', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
	    azelr = [azelr; az, el2, norm(s)];
        prn_act = [prn_act, prns(i)];
        % disp(prns(i))
    end
end

scatter3(r_gs(1)/1000, r_gs(2)/1000, r_gs(3)/1000, 'r');

%%
skyplot(azelr(:, 1)*180/pi, azelr(:, 2)*180/pi)

for i = 1:length(azelr)
    % skyplot(azelr(i, 1)*180/pi, azelr(i, 2)*180/pi);
    % hold on
    text(azelr(i, 1), azelr(i, 2)*180/pi, num2str(prn_act(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

%% 
figure
for i = 1:length(azelr)
    polarplot([0, azelr(i, 1)], [0, azelr(i, 2)], 'k')
    hold on
end

figure
for i = 1:length(azelr)
    polarplot([0, azelr(i, 1)], [0, pi/2 - azelr(i, 2)], 'k')
    text(azelr(i, 1), pi/2 - azelr(i, 2), num2str(prn_act(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    hold on
end

%% Values for image 
clc;
theta = azelr(:, 1)*180/pi;
r = pi/2 - azelr(:, 2);
r = r/max(r);

disp([theta, r*1500])