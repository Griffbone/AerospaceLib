% Function to calculate B-plane parameters of an object in a hyperbolic
% orbit 
% 
% Author: Griffin Jourda 10/14/22
% 
%	Inputs 
%		r		:	position 
%		v		:	velocity 
%		mu		:	gravitational parameter
% 
function [BT, BR, LTOF] = BPlane(r, v, mu) 
	eps = (norm(v)^2)/2 - mu/norm(r); 

	h = cross(r, v)/norm(cross(r, v));
	ev = (1/mu)*((norm(v)^2 - mu/norm(r))*r - dot(r, v)*v);
	
	e = norm(ev);
	a = -mu/(2*eps); 
	b = -a*sqrt(e^2 - 1); 
	del = acos(1/e); 
	
	S = (ev/e)*cos(del) + cross(h, ev)/norm(cross(h, ev))*sin(del);
	T = cross(S, [0 0 1])/norm(cross(S, [0 0 1]));
	R = cross(S, T); 
	B = b*cross(S, h);
	
	BT = dot(B, T); 
	BR = dot(B, R);
	
	LTOF = norm(B - r)/norm(v);
end
