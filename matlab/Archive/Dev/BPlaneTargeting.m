% Function to perform b-plane targeting in order to determine TCM maneuver
% requirements 
% 
% Author: Griffin Jourda 10/15/22
% 
%	Inputs 
%		jd_tcm	:	Julian date of the TCM maneuver
%		jd_a	:	estimated arrival Julian date (from solving Lambert's
%					problem)
%		r_tcm	:	heliocentric position at TCM maneuver (km/s)
%		v_tcm	:	heliocentric velocity at TCM maneuver (km/s)
%		BT_t	:	target B*T (km)
%		BR_t	:	target B*R (km)
%		pn		:	target planet (1 = Mercury, 2 = Venus, ... 8 = Neptune) 
%		
%	Outputs
%		dv		:	delta-v vector in heliocentric frame (km/s)

function [dv] = BPlaneTargeting(jd_tcm, jd_a, r_tcm, v_tcm, BT_t, BR_t, pn)
	% Get unperturbed B*T and B*R and dB/dV Jacobian
	[jac, BT1, BR1, ~] = BPlaneFDiff(r_tcm, v_tcm, jd_tcm, jd_a, pn); 
	dv = [0; 0; 0];
	n = 0;
	its = 4;

	% Perform corrections until B*T and B*R are within 1m of target
	while abs(BT_t - BT1) > 1e-3 && abs(BR_t - BR1) > 1e-3
		[BT1, BR1, ~] = propToBPlane(r_tcm, v_tcm + dv', jd_tcm, jd_a, pn);
		dv = dv + (jac^-1)*[BT_t - BT1; BR_t - BR1; 0];
		
		its = its + 1;
		n = n + 1;
	end

	% Targeting information/performance printout
	fprintf('B-Plane targeting converged in %i iterations and %i runs \n', n, its)
	fprintf('Final B*T: %.4f km \n', BT1)
	fprintf('Final B*R: %.4f km \n', BR1)
	fprintf('Final Delta-V: %.4f m/s \n', norm(dv)*1000)
end

function [dBdV, BT_u, BR_u, LTOF_u] = BPlaneFDiff(r, v, jd_tcm, jda, pn)
	dBdV = ones(3);

	% Propagate once without maneuver
	[BT1, BR1, LTOF1] = propToBPlane(r, v, jd_tcm, jda, pn);
	BT_u = BT1; 
	BR_u = BR1; 
	LTOF_u = LTOF1;

	% Vx perturbation
	[BT, BR, LTOF] = propToBPlane(r, v + [1/1000, 0, 0], jd_tcm, jda, pn);
	dBdV(1, 1) = (BT - BT1);
	dBdV(2, 1) = (BR - BR1);
	dBdV(3, 1) = (LTOF - LTOF1);

	% Vy perturbation
	[BT, BR, LTOF] = propToBPlane(r, v + [0, 1/1000, 0], jd_tcm, jda, pn);
	dBdV(1, 2) = (BT - BT1);
	dBdV(2, 2) = (BR - BR1);
	dBdV(3, 2) = (LTOF - LTOF1);

	% Vz perturbation
	[BT, BR, LTOF] = propToBPlane(r, v + [0, 0, 1/1000], jd_tcm, jda, pn);
	dBdV(1, 3) = (BT - BT1);
	dBdV(2, 3) = (BR - BR1);
	dBdV(3, 3) = (LTOF - LTOF1);
	
	dBdV = dBdV*1000;
end

function [BT, BR, LTOF] = propToBPlane(r0, v0, jd1, jd2, pn)
	load('constants.mat', 'mu_s');
	
	options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 2700, 'InitialStep', 60, 'Events', @(t, y) SOIEvent(t, y, jd1, pn));
	[t, y, ~, ~, ~] = ode45(@(t, y) twoBodyEOM(t, y, mu_s), [0, (jd2 - jd1)*86400], [r0', v0'], options);
	
	% Final satellite state 
	r = y(end, 1:3);
	v = y(end, 4:6); 
	
	% Final target planet state 
	[rp, vp, mu] = planetState(jd1 + max(t)/86400, pn);
	
	% B-plane parameters
	[BT, BR, LTOF] = BPlane(r - rp, v - vp, mu); 
end

