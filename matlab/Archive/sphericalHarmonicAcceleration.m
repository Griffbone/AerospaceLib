% Function to calculate gravity in the planet-fixed frame given latitude,
% longitude and spherical harmonic coefficients 
% 
% Author: Griffin Jourda 3/10/2023
% 
% Inputs 
%	phi	:	latitude (rad)
%	lam	:	longitude (rad)
% 	r	:	radius (m)
%	CS	:	Combined de-normalized C and S spherical harmonic coeffients
%	mu	:	central body gravitational parameter (m^3/s^2)
%	r_p	:	central body radius (m)
%
% Outputs
%	x	:	planet-fixed x-axis acceleration (m/s^2)
%	y	:	planet-fixed y-axis acceleration (m/s^2)
%	z	:	planet-fized z-axis acceleration (m/s^2)
%	V	:	gravitational potential at given point (m^2/s^2)

function [x, y, z, U] = sphericalHarmonicAcceleration(phi, lam, r, CS, mu, r_p)
	x = 0; 
	y = 0; 
	z = 0;
	
	mur2 = mu./(r_p.^2);
	nmax = size(CS, 1) - 1;
	U = 0;

	for n = 0:nmax
		inner = 0; 
		for m = 0:n
			Cnm = CS(n+1, m+1);

			if m == 0
				Snm = 0;
			else
				Snm = CS(m, n+1);
			end 

			inner = inner + lpoly(n, m, sin(phi)).*(Cnm.*cos(m.*lam) + Snm.*sin(m.*lam));

			if m == 0
				xnm = mur2.*(-Cnm.*V(n+1, 1, r_p, r, phi, lam));
				ynm = mur2.*(-Cnm.*W(n+1, 1, r_p, r, phi, lam));
			else
				Vnp1mp1 = V(n+1, m+1, r_p, r, phi, lam);
				Wnp1mp1 = W(n+1, m+1, r_p, r, phi, lam);
				Vnp1mm1 = V(n+1, m-1, r_p, r, phi, lam);
				Wnp1mm1 = W(n+1, m-1, r_p, r, phi, lam);

				xnm = mur2.*0.5.*((-Cnm.*Vnp1mp1 - Snm.*Wnp1mp1) + (factorial(n - m + 2)./factorial(n - m)).*(Cnm.*Vnp1mm1 + Snm.*Wnp1mm1));
				ynm = mur2.*0.5.*((-Cnm.*Wnp1mp1 + Snm.*Vnp1mp1) + (factorial(n - m + 2)./factorial(n - m)).*(-Cnm.*Wnp1mm1 + Snm.*Vnp1mm1));
			end

			znm = mur2.*(n - m + 1).*(-Cnm.*V(n+1, m, r_p, r, phi, lam) - Snm.*W(n+1, m, r_p, r, phi, lam));
			
			x = x + xnm;
			y = y + ynm;
			z = z + znm;
		end
		U = U + ((r_p/r).^n).*inner;
	end

	U = (mu./r).*U;
end

function [Vnm] = V(n, m, r_p, r, phi, lam)
	Vnm = ((r_p./r).^(n + 1)).*lpoly(n, m, sin(phi)).*cos(m.*lam);
end

function [Wnm] = W(n, m, r_p, r, phi, lam)
	Wnm = ((r_p./r).^(n + 1)).*lpoly(n, m, sin(phi)).*sin(m.*lam);
end