% Function to propagate an orbit by solving Kepler's equation. 
%
% Author: Griffin Jourda 10/01/22
% 
%	Inputs 
%		r	:	initial position 
%		v	:	initial velocity 
%		t	:	array of times to calculate satellite state
%		mu	:	gravitational parameter
% 
%	Outputs
%		t	:	array of satellite state was calculated at
% % % % % % % %		y_out	:	array of
%		rs	:	array of position vectors
%		vs	:	array of velocity vectors 
%		tas	:	array of true anomalies (rad)

function [t, rs, vs, tas] = keplerProp(r, v, t, mu)
	[a, e, i, raan, w, ta] = RV2COE(r, v, mu);
	
	rs = [r]; 
	vs = [v]; 
	tas = [ta]; 
	
	E = atan2(sin(ta)*sqrt(1 - e^2), e + cos(ta));
	M0 = E - e*sin(E);
	n = sqrt(mu/a^3);
	
	for ti = t(2:end)
		M = M0 + ti*n;
		[~, ta] = keplerE(M, e); 
		[r, v] = COE2RV(a, e, i, raan, w, ta, mu); 
		
		rs = [rs; r];
		vs = [vs; v];
		tas = [tas, ta];
	end
end