function [M] = cansatControl(t, y)
	M = (-1/10e3)*y(5:7);
end