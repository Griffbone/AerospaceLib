r0 = y(1, 1:3)/1000;
r0_gmat = y_gmat(1, 1:3);
disp(r0) 
disp(r0_gmat)

jd_utc = jd0;
[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp)';

r_ecef = R_icrs_itrf*r0';
disp(r_ecef')
disp(y_gmat_ecef(1, 1:3))


plot3(y_gmat_ecef(:, 1), y_gmat_ecef(:, 2), y_gmat_ecef(:, 3))

% r_ecef = []; 
% 
% for i = 1:size(y_gmat, 1) 
% 	ri = y_gmat(i, 1:3)';
% 
% 	jd_utc = jd0 + t_gmat(i)/86400; 
% 	[jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
% 	R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp)';
% 
% 	ri_ecef = R_icrs_itrf*ri;
% 
% 	r_ecef = [r_ecef; ri_ecef'];
% end
% 
% %%
% [t_gmat, y_gmat_ecef, jd0] = interpGMAT('ephem/sso_twoBodySHGrav_ecef.e');
% 
% plot(t_gmat, y_gmat_ecef(:, 1:3) - r_ecef(:, 1:3))
% hold on 
% yline(6378)
% 
% %% 
% % jd_utc = greg2jd(2023, 1, 1, 0, 0, 0)
% % EOPData = loadEOPData('EOPdata.txt');
% % [jd_ut1, ~, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPData);
% % 
% % R_icrs_itrf = eci2ecef_iau(jd_tt, jd_ut1, xp, yp)'