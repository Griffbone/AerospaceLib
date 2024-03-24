% Function to calculate local body-fixed gravity using a spherical
% harmonics model 
% 
% Author: Griffin Jourda 4/01/2023
%
% Inputs
%	r	:	body-centered position vector
%	C	:	denormalized Cnm coefficients 
%	S	:	denormalized Snm coefficients 
%	rp	:	planet radius 
%	mu	:	planet gravitational parameter
% 
% Outputs 
%	xdd		:	body-fixed x-axis acceleration 
%	ydd		:	body-fixed y-axis acceleration 
%	zdd		:	body-fixed z-axis acceleration

function [xdd, ydd, zdd] = sphericalHarmonicGravity(r, C, S, rp, mu)
	nmax = size(C, 1) - 1;
	[V, W] = getVW(r, rp, nmax);
	nmax = nmax - 1;

	murp2 = mu/(rp*rp);
	murp22 = murp2/2; 
	
	xdd = 0; 
	ydd = 0;
	zdd = 0;
	for n = 0:nmax
		nidx = n + 1;
		for m = 0:n
			midx = m + 1; 
			
			if m == 0
				xdd = xdd + (murp2)*(-C(nidx, 1)*V(nidx+1, 2));
				ydd = ydd + (murp2)*(-C(nidx, 1)*W(nidx+1, 2));
			else
				xdd = xdd + (murp22)*( ...
					(-C(nidx, midx)*V(nidx+1, midx+1) - S(nidx, midx)*W(nidx+1, midx+1)) + ...
					(factorial(n - m + 2)/factorial(n - m))*(C(nidx, midx)*V(nidx+1, midx-1) + ...
					S(nidx, midx)*W(nidx+1, midx-1)));
				
				ydd = ydd + (murp22)*( ...
					(-C(nidx, midx)*W(nidx+1, midx+1) - S(nidx, midx)*V(nidx+1, midx+1)) + ...
					(factorial(n - m + 2)/factorial(n - m))*(C(nidx, midx)*W(nidx+1, midx-1) + ...
					S(nidx, midx)*V(nidx+1, midx-1)));
			end
			
			
			zdd = zdd + (murp2)*((n - m + 1)*(-C(nidx, midx)*V(nidx + 1, midx) - S(nidx, midx)*W(nidx + 1, midx + 1)));
		end
	end
	
end

function [V, W] = getVW(r, rp, nmax)
	x = r(1); y = r(2); z = r(3); 
	r2 = x*x + y*y + z*z;
	rp2 = rp*rp;

	zrpr2 = z*rp/r2;
	rp2r2 = rp2/r2;
	
	V = zeros(nmax+1);
	W = zeros(nmax+1);
	
	for m = 0:nmax
		midx = m+1;
	
		[Vmm, Wmm] = vwmm(m, r, rp);
		V(midx, midx) = Vmm; 
		W(midx, midx) = Wmm;
	
		for n = m+1:nmax
			nidx = n+1;
	
			if n - m < 2
				Vprv = 0;
				Wprv = 0;
			else 
				Vprv = V(nidx-2, midx);
				Wprv = W(nidx-2, midx);
			end
			
			a = ((2*n - 1)/(n - m))*(zrpr2);
			b = ((n + m - 1)/(n - m))*(rp2r2);

			V(nidx, midx) = a*V(nidx - 1, midx) - b*Vprv;
			W(nidx, midx) = a*W(nidx - 1, midx) - b*Wprv;
		end
	end
end

function [Vmm, Wmm] = vwmm(m, r, rp)
	x = r(1); y = r(2); z = r(3); 
	r2 = (x*x +y*y + z*z);
	xrpr2 = x*rp/r2;
	yrpr2 = y*rp/r2;

	if m == 0
		Vmm = rp/norm(r);
		Wmm = 0; 
	else
		[Vm1, Wm1] = vwmm(m-1, r, rp);

		Vmm = (2*m - 1)*(xrpr2*Vm1 - yrpr2*Wm1);
		Wmm = (2*m - 1)*(xrpr2*Wm1 + yrpr2*Vm1);
	end
end