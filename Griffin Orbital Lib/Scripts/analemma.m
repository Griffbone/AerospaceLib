clc; clear;
EOPData = loadEOPData('EOPdata.txt');

%% 
jd0 = greg2jd(2022, 1, 1, 12, 0, 0);
jdf = greg2jd(2023, 1, 1, 12, 0, 0);

azel = []; 

lon = 0.0015; 
lat = 51.48;
tspan = jd0:jdf;

for jd = tspan
	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd, EOPData); 
	R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp);

	s = sunVector(jd);
	s_ecef = R_icrs_itrf*s;

	s_enu = ECEF2ENU_simple(s_ecef, lat, lon);
	[az, el] = ENU2azel(s_enu);

	azel = [azel; az, el];
end

%% Plotting
% skyplot(azel(:, 1), azel(:, 2))
scatter(azel(:, 1), azel(:, 2), 10, tspan, 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2)
xlabel('Azimuth [deg]') 
ylabel('Elevation [deg]')
colorbar()
grid on


i_ss = find(azel(:, 2) == max(azel(:, 2)));
jd_ss = tspan(i_ss); 
[~, ~, ~, ~, ~, ~, ss_date] = jd2greg(jd_ss);

i_ws = find(azel(:, 2) == min(azel(:, 2)));
jd_ws = tspan(i_ws); 
[~, ~, ~, ~, ~, ~, ws_date] = jd2greg(jd_ws);

hold on

scatter(azel(i_ss, 1), azel(i_ss, 2), 'k')
scatter(azel(i_ws, 1), azel(i_ws, 2), 'k', 'square')
legend({'-', sprintf('Summer Solstice: %s', ss_date(1:11)), sprintf('Winter Solstice: %s', ws_date(1:11))})

phi = min(azel(:, 2)) + (max(azel(:, 2)) - min(azel(:, 2)))/2;
[~, order] = sort(abs(azel(:, 2) - phi));

i_s1 = order(1); 
i_s2 = order(2);

scatter(azel(i_s1, 1), azel(i_s1, 2), 'k')
scatter(azel(i_s2, 1), azel(i_s2, 2), 'k')
yline(phi)

jd_s1 = tspan(i_s1); 
[~, ~, ~, ~, ~, ~, s1_date] = jd2greg(jd_s1);

jd_s2 = tspan(i_s2); 
[~, ~, ~, ~, ~, ~, s2_date] = jd2greg(jd_s2);

disp(s1_date) 
disp(s2_date)








