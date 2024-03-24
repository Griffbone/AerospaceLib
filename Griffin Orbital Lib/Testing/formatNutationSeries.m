clc; clear; 

f = fopen('nutationSeries.txt');
f2 = fopen('formattedSeries.txt', 'w');

line = fgetl(f);
while ischar(line)
	line(line == ' ') = ',';
	line = ['[', line(1:end-1), '],\n'];

	fprintf(f2, line);

	line = fgetl(f);
end

fclose(f);
fclose(f2);