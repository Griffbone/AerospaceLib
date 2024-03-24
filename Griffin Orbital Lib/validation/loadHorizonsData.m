function [data] = loadHorizonsData(fname)
	fid = fopen(fname);
	line = fgetl(fid);
	data = [];

	while ischar(line)
		ldat = split(line, ',');
		
		data = [data; str2double(ldat{1}), str2double(ldat{3}), str2double(ldat{4}), str2double(ldat{5}), str2double(ldat{6})];
	
		line = fgetl(fid);
	end
	fclose(fid);

	data(:, 1) = data(:, 1) - data(:, 2)/86400;

	data = [data(:, 1), data(:, 3:end)*1000];
end