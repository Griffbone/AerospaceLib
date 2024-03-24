% Solve Kepler's problem for eccentric anomaly
%
% Author: Griffin Jourda 3/25/2022
%
%	Inputs
%		M		:	Mean anomaly (rad)
%		e		:	eccentricity
%
%	Ouputs
%		E		:	eccentric anomaly (rad)
%		ta		:	true anomaly (rad)

function [E, ta] = keplerE(M, e)
	if M < pi
		E = M + e/2; 
	else
		E = M - e/2;
	end

	diff = 1;
	
	while abs(diff) > 1e-9
		E0 = E;
		dMdE = 1 - e*cos(E0);
		E = E0 - (E0 - e*sin(E0) - M)/dMdE;
		diff = E - E0;
	end	

	ta = wrapTo2Pi(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2)));
end
