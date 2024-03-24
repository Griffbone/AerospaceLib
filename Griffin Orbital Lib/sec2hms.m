function [h, m, s] = sec2hms(sec) 
	h = floor(sec/3600); 
	m = floor((sec - h*3600)/60);
	s = (sec - h*3600 - m*60);
end