function [r] = interpHorizonsData(data, jd_utc)
	r = interp1(data(:, 1), data(:, 2:end), jd_utc, 'spline')';
end