% Convert orbital elements to inertial position and velocity
%
% Author: Griffin Jourda 3/30/2022
%
%	Inputs
%		a		:	semimajor axis
%		e		:	eccentricity 
%		i		:	inclination (rad) 
%		raan	:	right ascension of the ascending node (rad) 
%		w		:	argument of periapsis (rad)
%		ta		:	true anomaly (rad) 
%		mu		:	gravitational parameter
% 		
% 	Outputs
% 		rvec	:	ECI position vector 
%		vvec	:	ECI velocity vector

function [rvec, vvec] = COE2RV(a, e, i, raan, w, ta, mu)
	P = a*(1 - e^2); 
	
	r_pqw = [P*cos(ta)/(1 + e*cos(ta)); P*sin(ta)/(1 + e*cos(ta)); 0];
	v_pqw = [-sqrt(mu/P)*sin(ta); sqrt(mu/P)*(e + cos(ta)); 0];
	
	rvec = R3(-raan)*R1(-i)*R3(-w)*r_pqw;
	vvec = R3(-raan)*R1(-i)*R3(-w)*v_pqw;

	rvec = rvec';
	vvec = vvec';
end