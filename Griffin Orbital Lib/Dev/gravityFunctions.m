clc; clear; 
[C, S] = loadSHData('EGM2008.txt', 21);


%% 
rp = 0.63781363e7;
mu = 0.3986004415e15;

n = 100;
phi = linspace(pi/2, -pi/2, n);
lam = linspace(0, 2*pi, n);
g = zeros(length(phi));
g_j2 = zeros(length(phi));

for i = 1:length(phi)
	lat = phi(i);
	
	for j = 1:length(lam)
		lon = lam(j);
		r = (6378e3 + 100e3)*[cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];
		
		[x, y, z] = sphericalHarmonicGravity(r, C, S, rp, mu);
		
		g(i, j) = dot([x; y; z], r/norm(r));
	end	
end

%% 
pcolor(g)
shading('interp')
colorbar

%% 
clc;
n = 200;

r = rand(n, 3);
r = r./vecnorm(r, 2, 2);
r = r*6378e3 + 500e3 + rand(n, 1)*500e3;

rp = 0.63781363e7;
mu = 0.3986004415e15;

e = zeros(n, 3);
for i = 1:n
	ri = r(i, :)';
	[x, y, z] = sphericalHarmonicGravity(ri, C, S, rp, mu);
	[xx, yy, zz] = gravitysphericalharmonic(ri', 20);
	
	e(i, :) = [abs(xx-x)*100000, abs(yy-y)*100000, abs(zz-z)*100000];
end

%% 
alt = (vecnorm(r, 2, 2) - 6378e3)/1000;

figure
histogram(alt, 20)
xlabel('Altitude [km]')
ylabel('Number of Occurance [-]')

figure; hold on;
histogram(e(:, 1), 20)
histogram(e(:, 2), 20)
histogram(e(:, 3), 20)
xlabel('Error [mgal]')
ylabel('Number of Occurances [-]')
legend('X', 'Y', 'Z')

figure; subplot(1, 3, 1);
scatter(alt, e(:, 1))

subplot(1, 3, 2)
scatter(alt, e(:, 2))

subplot(1, 3, 3)
scatter(alt, e(:, 2))

%% Functions
function [xdd, ydd, zdd] = sphericalHarmonicGravity(r, C, S, rp, mu)
	nmax = size(C, 1) - 1;
	[V, W] = getVW(r, rp, nmax);
	
	nmax = size(C, 1) - 2;
	
	xdd = 0; 
	ydd = 0;
	zdd = 0;
	for n = 0:nmax
		nidx = n + 1;
		for m = 0:n
			midx = m + 1; 
			
			if m == 0
				xdd = xdd + (mu/rp^2)*(-C(nidx, 1)*V(nidx+1, 2));
				ydd = ydd + (mu/rp^2)*(-C(nidx, 1)*W(nidx+1, 2));
			else
				xdd = xdd + (mu/(2*rp^2))*( ...
					(-C(nidx, midx)*V(nidx+1, midx+1) - S(nidx, midx)*W(nidx+1, midx+1)) + ...
					(factorial(n - m + 2)/factorial(n - m))*(C(nidx, midx)*V(nidx+1, midx-1) + ...
					S(nidx, midx)*W(nidx+1, midx-1)));
				
				ydd = ydd + (mu/(2*rp^2))*( ...
					(-C(nidx, midx)*W(nidx+1, midx+1) - S(nidx, midx)*V(nidx+1, midx+1)) + ...
					(factorial(n - m + 2)/factorial(n - m))*(C(nidx, midx)*W(nidx+1, midx-1) + ...
					S(nidx, midx)*V(nidx+1, midx-1)));
			end
			
			
			zdd = zdd + (mu/rp^2)*((n - m + 1)*(-C(nidx, midx)*V(nidx + 1, midx) - S(nidx, midx)*W(nidx + 1, midx + 1)));
		end
	end
	
end

function [V, W] = getVW(r, rp, nmax)
	x = r(1); y = r(2); z = r(3); 
	r2 = x*x + y*y + z*z;
	rp2 = rp*rp;
	
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
			
			a = ((2*n - 1)/(n - m))*(z*rp/r2);
			b = ((n + m - 1)/(n - m))*(rp2/r2);

			V(nidx, midx) = a*V(nidx - 1, midx) - b*Vprv;
			W(nidx, midx) = a*W(nidx - 1, midx) - b*Wprv;
		end
	end
end

function [Vmm, Wmm] = vwmm(m, r, rp)
	x = r(1); y = r(2); z = r(3); 
	r2 = (x*x +y*y + z*z);

	if m == 0
		Vmm = rp/norm(r);
		Wmm = 0; 
	else
		[Vm1, Wm1] = vwmm(m-1, r, rp);

		Vmm = (2*m - 1)*((x*rp/r2)*Vm1 - (y*rp/r2)*Wm1);
		Wmm = (2*m - 1)*((x*rp/r2)*Wm1 + (y*rp/r2)*Vm1);
	end
end