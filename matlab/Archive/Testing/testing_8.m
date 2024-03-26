clc; clear; 
addpath('C:\Users\Griffin\OneDrive\GaTech\Griffin Orbital Lib')

fname = 'eopData.txt';
[jd, xp, yp, dut1, dtai] = loadEOPData(fname);
EOPDat = [jd, xp, yp, dut1, dtai];

%% Testing
jd_utc = greg2jd(1999, 3, 4, 0, 0, 0) - 13/86400;
[jd_ut1, jd_tai, jd_tt, xp, yp] = interpEOPData(jd_utc, EOPDat);

R_icrs_itrs = ICRS2ITRS(jd_tt, jd_ut1, xp, yp);
R_itrs_icrs = R_icrs_itrs';

r_wgs = [19440.953805, 16881.609273, -6777.115092]';

r_icrs = R_itrs_icrs*r_wgs