% Function to evaluate the Legendre polynomial of order n and degree m for
% an input of t
% 
% Author: Griffin Jourda 3/10/2023
% 
% Inputs 
%	n	:	degree (1-inf)
%	m	:	order (0-n)
%	t	:	value to evaluate polynomial at
%
% Outputs
%	Pnm	:	n, m Legendre polynomial evaluted at t

function [Pnm] = lpoly(n, m, t)
	sigma = 0;
	r = floor((n - m)./2);

	parfor k = 0:r
		sigma = sigma + ...
		((-1).^k).*factorial(2.*n - 2.*k)./(factorial(k).*factorial(n - k).*factorial(n - m - 2.*k)).*t.^(n - m - 2*k);
	end

	Pnm = ((2.^(-n)).*(1 - t.^2).^(m./2)).*sigma;
end