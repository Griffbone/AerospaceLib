% Function to load spherical harmonics coefficients from a data file
% 
% Author: Griffin Jourda 3/28/2023
% 
% Inputs
%	fname	:	file name of gravity model coefficients 
%	maxdeg	:	maximum degree of coefficients to return 
% Outputs 
%	C		:	array of denormalized C coefficients 
%	S		:	array of denormalized S coefficients 

function [C, S] = loadSHData(fname, maxdeg)
	data = readmatrix(fname, 'FileType', 'text');
	n = floor(data(:, 2)); 
	m = floor(data(:, 3));
	C_norm = data(:, 4);
	S_norm = data(:, 5);

	C = zeros(maxdeg + 1);
	S = zeros(maxdeg + 1); 

	for i = 1:length(n) 
		ni = n(i);
		mi = m(i);
		Cbar = C_norm(i);
		Sbar = S_norm(i); 

		if mi == 0
			k = 1; 
		else 
			k = 2; 
		end
		
		C(ni + 1, mi + 1) = Cbar/sqrt(factorial(ni + mi)/(k*(2*ni + 1)*factorial(ni - mi)));
		S(ni + 1, mi + 1) = Sbar/sqrt(factorial(ni + mi)/(k*(2*ni + 1)*factorial(ni - mi)));
		
		if (ni == maxdeg) && (mi == ni)
			break 
		end
	end
end