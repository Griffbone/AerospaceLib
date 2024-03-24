% Function to calculate gravity anomaly in the planet-fixed coordinate
% frame given latitude, longitude, and spherical harmonic coefficients.
% Gravity anomaly is defined as the gravitational acceleration minus the
% contributions from C00 and C20.
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

function [x, y, z, x_lt, y_lt, z_lt] = sphericalHarmonicAnomaly(phi, lam, r, CS, mu, r_p)
	x = 0; 
	y = 0; 
	z = 0;

	x_lt = 0;
	y_lt = 0; 
	z_lt = 0;
	
	mur2 = mu./r_p.^2;
	nmax = size(CS, 1) - 1;

	for n = 0:nmax
		for m = 0:n
			Cnm = CS(n+1, m+1);

			if m == 0
				Snm = 0;
			else 
				Snm = CS(m, n+1);
			end 

			if m == 0
				xnm = mur2.*(-Cnm.*V(n+1, 1, r_p, r, phi, lam));
				ynm = mur2.*(-Cnm.*W(n+1, 1, r_p, r, phi, lam));
			else
				Vnp1mp1 = V(n, m, r_p, r, phi, lam);
				Wnp1mp1 = W(n, m, r_p, r, phi, lam);
				Vnp1mm1 = V(n+1, m-1, r_p, r, phi, lam);
				Wnp1mm1 = W(n+1, m-1, r_p, r, phi, lam);

				xnm = mur2.*0.5.*((-Cnm.*Vnp1mp1 - Snm.*Wnp1mp1) + (factorial(n - m + 2)./factorial(n - m)).*(Cnm.*Vnp1mm1 + Snm.*Wnp1mm1));
				ynm = mur2.*0.5.*((-Cnm.*Wnp1mp1 + Snm.*Vnp1mp1) + (factorial(n - m + 2)./factorial(n - m)).*(-Cnm.*Wnp1mm1 + Snm.*Vnp1mm1));
			end

			znm = mur2.*((n - m+1).*(-Cnm.*V(n+1, m, r_p, r, phi, lam) - Snm.*W(n+1, m, r_p, r, phi, lam)));
			
			if (n == 0 && m == 0) || (n == 2 && m == 0)
				x_lt = x_lt + xnm;
				y_lt = y_lt + ynm;
				z_lt = z_lt + znm;
			end
			
			x = x + xnm;
			y = y + ynm;
			z = z + znm;
		end
	end
end

function [Vnm] = V(n, m, r_p, r, phi, lam)
	Vnm = ((r_p./r).^(n + 1)).*lpoly(n, m, sin(phi)).*cos(m.*lam);
end

function [Wnm] = W(n, m, r_p, r, phi, lam)
	Wnm = ((r_p./r).^(n + 1)).*lpoly(n, m, sin(phi)).*sin(m.*lam);
end