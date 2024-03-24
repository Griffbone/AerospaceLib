% Function to convert an angle in degrees to hours, minutes, and seconds 
% 
% Author: Griffin Jourda 10/06/22
% 
%	Inputs
%		deg	:	angle in degrees 
%	
%	Outputs 
%		h	:	hour angle 
%		m	:	minute angle 
%		s	:	second angle 

function [h, m, s] = deg2hms(deg)
	hrs = deg/15; 

	h = floor(hrs); 
	m = floor((hrs - h)*60);
	s = (hrs - h - m/60)*3600;
end