clc; clear;

addpath('C:\Users\Griffin\OneDrive\GaTech\Griffin Orbital Lib')

jd_utc = greg2jd(2018.6, 2, 1, 0, 0, 0);
EOPDATA = loadEOPData('eopData.txt');
[jd_ut1, jd_tai, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPDATA);

r_eci = [6378; 0; 0];
r_ecef = ECI2ECEF_simple(r_eci, jd_ut1);
r_itrf = ICRS2ITRS(jd_tt, jd_ut1, xp, yp)*r_eci;

disp(norm(r_ecef - r_itrf)) 
disp(acosd(dot(r_ecef, r_itrf)/(norm(r_ecef)*norm(r_itrf))))