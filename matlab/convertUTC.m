% Function to convert UTC time to various other time systems, given year,
% month, day, and UTC H:M:S. Returns UT1, TAI, TT, and TDB times in seconds
% corresponding to given UTC HMS (can be converted to hms using sec2hms
% function). 
% 
% Author: Griffin Jourda 12/23/2022
% 
%	Inputs 
%		yr	:	UTC year 
%		mo	:	UTC month 
%		day	:	UTC day
%		h	:	UTC hour 
%		m	:	UTC minute
%		s	:	UTC second 
%		DAT	:	difference between UTC and TAI (accumulated leap seconds)
%		DUT1:	difference between UTC and UT1
% 
%	Outputs
%		UT1	:	UT1 seconds corresponding to UTC hms
%		TAI	:	atomic time seconds corresponding to UTC hms
%		TT	:	terrestrial time seconds corresponding to UTC hms
%		TDB :	barycentric dynamical time corresponding to UTC hms
function [UT1, TAI, TT, TDB] = convertUTC(yr, mo, day, h, m, s, DUT1, DAT)
	% Convert UTC time to seconds
	UTC = h*3600 + m*60 + s;

	% Find UT1, TAI, and TT
	UT1 = UTC + DUT1;
	TAI = UTC + DAT; 
	TT = TAI + 32.184; 
	
	% Find JD in TT and TT Julian Centuries
	[h_TT, m_TT, s_TT] = sec2hms(TT); 
	jd_TT = greg2jd(yr, mo, day, h_TT, m_TT, s_TT);
	T = (jd_TT - 2451545)/36525;

	% Convert TT to TDB ~ Teph (Kaplan, 2005 eq 2.6)
	TDB = TT + 0.001657*sin(628.3076*T + 6.2401) + ...
		0.000022*sin(575.3385*T + 4.2970) + ...
		0.000014*sin(1256.6152*T + 6.1969) + ...
		0.000005*sin(606.9777*T + 4.0212) + ...
		0.000005*sin(52.9691*T + 0.4444) + ...
		0.000002*sin(21.3299*T + 5.5431) + ...
		0.000010*T*sin(628.3076*T + 4.2490);
end