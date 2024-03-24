% Convert a Gregorian date to its corresponding Julian date 
%
% Author: Griffin Jourda 2/13/2022
%
%	Inputs
%		year	:	gregorian calendar year
%		month	:	gregorian calendar month
%		day		:	gregorian calendar day
%		hr		:	UT hour
%		min		:	UT minute	
%		sec		:	UT second
%
%	Ouputs
%		jdn		:	Julian day number
%		jd_frac :	Julian day fraction (days)

function [jd] = greg2jd(year, month, day, hr, min, sec)	
	jdn = 367*year - floor((7/4) * (year + floor((month + 9)/12))) + floor((275/9)*month) + day + 1721013.5;
	jd_frac = (sec/60/60 + min/60 + hr)/24;

	jd = jdn + jd_frac;
end