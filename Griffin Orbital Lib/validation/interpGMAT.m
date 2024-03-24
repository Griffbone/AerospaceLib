function [t, y, epoch] = interpGMAT(fname)
	fid = fopen(fname);
	line = fgetl(fid); 
	header = 1;

	data = [];

	while ischar(line) 

		if contains(line, 'ScenarioEpoch')
			[~, epoch] = strtok(line, ' ');	
			epoch = strtrim(epoch);
			epoch = split(epoch, ' ');

			d = str2double(epoch{1});
			m = month2num(epoch{2});
			y = str2double(epoch{3});
			hms = epoch{4};

			h = str2double(hms(1:2));
			mi = str2double(hms(4:5));
			s = str2double(hms(7:end));

			epoch = greg2jd(y, m, d, h, mi, s);
		end

		if contains(line, 'EphemerisTimePosVel')
			header = 0;
		end

		if ~header && ~contains(line, 'EphemerisTimePosVel')
			c = strsplit(line);
			
			if length(c) == 7
				datai = [str2double(c{1}), str2double(c{2}), str2double(c{3}), ...
					str2double(c{4}), str2double(c{5}), str2double(c{6}), str2double(c{7})];
				data = [data; datai];
			end
		end

		line = fgetl(fid); 
	end
	fclose(fid);

	t = data(:, 1);
	y = data(:, 2:end);
end

function [m] = month2num(month)
	switch lower(month) 
		case 'jan'
			m = 1;
		case 'feb'
			m = 2;
		case 'mar'
			m = 3;
		case 'apr'
			m = 4;
		case 'may'
			m = 5;
		case 'jun'
			m = 6;
		case 'jul'
			m = 7;
		case 'aug' 
			m = 8; 
		case 'sep'
			m = 9;
		case 'oct'
			m = 10;
		case 'nov'
			m = 11; 
		case 'dec'
			m = 12;
	end
end