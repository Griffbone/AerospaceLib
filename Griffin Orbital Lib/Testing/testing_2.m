clc; clear;

jd = greg2jd(2022, 10, 17, 4, 0, 0);
[rp, vp, ~] = planetState(jd, 3);

r0 = rp + 6378*rp/norm(rp); 
v0 = vp + 8*vp/norm(vp);

t = linspace(0, 86400*365, 1000);
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
[t, y] = ode45(@(t, y) nBodyEOM(t, y, jd, [3]), t, [r0'; v0'], options);


plotTransfer(jd, t, y, [1, 2, 3, 4])


function [ydot] = nBodyEOM(t, y, j0, ps)
	jd = j0 + t/86400;
	
	r = y(1:3); 
	v = y(4:6); 
	a = -(132712440018./(norm(r)^3)).*r;
	
	for i = 1:length(ps)
		pn = ps(i); 
		[rp, ~, mu] = planetState(jd, pn);
		rsp = r - rp'; 
		rspn = norm(rsp);
		a = a - (mu./rspn^3).*rsp;
		disp(t/86400)
	end
	
	ydot = [v; a];
end