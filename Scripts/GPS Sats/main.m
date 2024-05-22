clc; clear; 
addpath('C:\Users\Griffin\Documents\AerospaceLib\matlab');

%% Propagate each TLE to date 
jd0 = greg2jd(2024, 5, 22, 0, 0, 0);
jdf = greg2jd(2024, 7, 8, 23, 0, 0);

r_ends = [];
prns = [];

for i = 1:31
	fid = fopen(['tles\tle_', num2str(i), '.txt']);
	l0 = fgetl(fid);
	l1 = fgetl(fid);
	l2 = fgetl(fid);
	fclose(fid);

	prns = [prns, str2num(regexp(l0, '(?<=PRN\s)\d{2}', 'match', 'once'))];
	disp(prns(i))
	disp(l0)

	[~, y_eci, y_ecef, ~, ~] = propTLE(l1, l2, jd0, jdf, 'J2');
	r_ends = [r_ends; y_ecef(end, 1:3)/1000];
end 

%% 
r_gs_ecef = LLA2ECEF(35.347566370233956, -117.8094102639529, 1068.0192)/1000;
plotInertialAxes(); 

for i = 1:31
	r = r_ends(i, :);
	s = r - r_gs_ecef;
	
	el = pi/2 - acos(dot(r_gs_ecef/norm(r_gs_ecef), s/norm(s)));

	if el > 0
		plot3([r_gs_ecef(1), r(1)], [r_gs_ecef(2), r(2)], [r_gs_ecef(3), r(3)], 'k')
		scatter3(r(1), r(2), r(3), 'k', 'filled', 'square');
		text(r(1), r(2), r(3), num2str(prns(i)), 'color', 'red', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
	end
end

scatter3(r_gs_ecef(1), r_gs_ecef(2), r_gs_ecef(3), 'r');
