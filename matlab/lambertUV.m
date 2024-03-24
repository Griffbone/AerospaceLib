% Function to solve Lambert's problem with universal variables
% 
% Author: Griffin Jourda 10/12/22
% 
%	Inputs 
%		r1	:	first position vector 
%		r2	:	second position vector 
%		tof :	time of flight 
%		mu	:	gravitational parameter 
%		tm	:	transfer method (1 = long, 2 = short)
% 
%	Outputs 
%		v1	:	initial velocity 
%		v2	:	final velocity 

function [v1, v2] = lambertUV(r1, r2, tof, mu, tm)
	flag = 0; 
	
	r1n = norm(r1); 
	r2n = norm(r2); 
	
	cdt = dot(r1, r2)/(r1n*r2n); 
	A = tm*sqrt(r1n*r2n*(1 + cdt)); 
	
	zlow = -4*pi; 
	zhigh = 4*pi^2;
	z = (zlow + zhigh)/2; 
	dt = 1;
	n = 0; 
	
	while abs(dt) > 1e-6
		[c, s] = stumpff(z); 
		yn = r1n + r2n + A*(z*s - 1)/sqrt(c); 
		
		if (A > 0) && (yn < 0)
			while yn < 0
				zlow = zlow + 0.1; 
				z = (zlow + zhigh)/2; 
				[c, s] = stumpff(z); 
				yn = r1n + r2n + A*(z*s - 1)/sqrt(c);
			end
		end 
		
		xn = sqrt(yn/c); 
		t = (1/sqrt(mu))*((xn^3)*s + A*sqrt(yn)); 
		
		if t < tof
			zlow = z; 
		else
			zhigh = z; 
		end
		
		z = (zlow + zhigh)/2; 
		dt = t - tof; 
		n = n + 1; 		
		
		if n > 50 
			flag = 1;
			break 
		end
	end
	
	if flag == 1
		v1 = [nan, nan, nan]; 
		v2 = [nan, nan, nan];
	else
		[c, s] = stumpff(z); 
		yn = r1n + r2n + A*(z*s - 1)/sqrt(c); 

		f = 1 - yn/r1n; 
		g = A*sqrt(yn/mu);
		gdot = 1 - yn/r2n; 

		v1 = (r2 - f*r1)/g;
		v2 = (gdot*r2 - r1)/g; 

		h1 = cross(r1, v1); 

		% Run again if solution is retrograde
		if h1(3) < 0
			[v1, v2] = lambertUV(r1, r2, tof, mu, -tm);
		end
	end
end