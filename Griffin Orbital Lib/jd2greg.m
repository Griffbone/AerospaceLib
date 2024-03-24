% Function to determine Gregordian calendar date from Julian date 
% 
% Author: Griffin Jourda 9/27/22
% 
%	Inputs
%		jd	:	Julian date 
%
%	Outputs
%		y	:	year
%		m	:	month
%		d	:	day
%		h	:	hour 
%		m	:	minute 
%		s	:	s
%		str	:	yyy-mm-dd hh:mm:ss.ssss string 

function [y, m, d, h, minute, s, str] = jd2greg(jd)
	lmonth = [31 28 31 30 31 30 31 31 30 31 30 31];

	% Determine year and if it is a leap year
	T1900 = (jd - 2415019.5)/365.25;
	y = 1900 + floor(T1900);

	leapyrs = floor((y - 1900 - 1)*0.25);
	days = (jd - 2415019.5) - ((y - 1900)*365 + leapyrs);

	if days < 1 
		y = y - 1;
		leapyrs = floor((y - 1900 - 1)*0.25); 
		days = (jd - 2415019.5) - ((y - 1900)*365 + leapyrs);
	end 

	if mod(y, 4) == 0
		lmonth(2) = 29;
	end

	dayofyr = floor(days);
	
	% Determine month and day (basically taken straight from Vallado's code
	% - I couldn't get this working for forever)
	m = 1; 
	sigma = 0; 	
	while sigma + lmonth(m) < dayofyr
		sigma = sigma + lmonth(m); 
		m = m + 1; 
	end
	
	d = dayofyr - sigma;
	
	% Determine hours, minutes, and seconds
	tau = (days - dayofyr)*24;
	h = floor(tau); 

	minute = floor((tau - h)*60);
	h = floor(tau);
	s = (tau - h - minute/60)*3600;

	% Make a date string 
	mons = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
	str = sprintf('%.2i %s %.2i %.2i:%.2i:%.3f UTC \n', d, mons{m}, y, h, minute, s);				% GMAT format
% 	str = sprintf('%i-%.2i-%.2i %.2i:%.2i:%.4f UTC \n', y, m, d, h, minute, s);						% Sane format
end