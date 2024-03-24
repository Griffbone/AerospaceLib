clc; clear; 
addpath('C:\Users\Griffin\OneDrive\GaTech\AE 4361\Library')
load('Constants.mat')

I = [28.4613, -3.2971, -2.7267;
	-3.2971, 17.6217, 1.8887;
	-2.7267, 1.8887, 13.9170];

r0 = [6500e3; 0; 0];				% 1-3
v0 = [0; 4500; 8000];				% 4-6
w0 = deg2rad([1; -7; 5]);			% 7-9
q0 = v2q(deg2rad([15; 30; 45]));	% 10-13
y0 = [r0; v0; w0; q0];
jd0 = greg2jd(2022, 11, 23, 0, 0, 0); 

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
tic()
[t, y] = ode45(@(t, y) fullStateEOM(t, y, I, jd0), linspace(0, 10e3, 10e3), y0, options);
toc()

%% Functions
% Functiont to calculate full 6DOF equations of motion for a satellite in
% LEO. Uses quaternion attitude representation and assumes torque-free
% motion. 
% 
% Author: Griffin Jourda 11/23/22
% 
%	Inputs 
%		t	:	time (s) 
%		y	:	satellite state [r; v; w; q] (m, m/s, rad/s, --) 
%		I	:	moment of inertia tensor (kg*m^2)
%		jd0	:	Julian date of propagation start
%
%	Outputs 
%		ydot	:	state derivative [v; a; wdot; qdot] (m/s, m/s^2,
%					rad/s^2, --)
function [ydot] = fullStateEOM(t, y, I, jd0)
	r = y(1:3); 
	v = y(4:6); 
	w = y(7:9); 
	q = y(10:13); 

	r_ecef = ECI2ECEF_simple(r, jd0 + t/86400);
	[lat, lon, alt] = ECEF2LLA(r_ecef, 'gd');


% 	B_ned = igrfmagm(alt, lat, lon, 2022.5)';			% date should not be hard coded but I'm too lazy to fix this RN
	B_ned = [30e-9; 0; 0];
	B_ecef = NED2ECEF(B_ned, lat, lon); 
	B_eci = R3(deg2rad(jd2gmst(jd0 + t/86400)))'*B_ecef;
	B_bf = q2R(q)*B_eci;

	% Orbital Dynamics 
	a = -(398600435507000/norm(r)^3)*r;

	% Attitude Dynamics
	wdot = I\(-cross(w, I*w));
	qdot = (1/2)*qMult([w; 0], q); 

	% Create state vector derivative
	ydot = [v; a; wdot; qdot];
end

function [v_ecef] = NED2ECEF(v_ned, lat, lon) 
	lat = deg2rad(lat); 
	lon = deg2rad(lon); 

	v_enu = [0 1 0; 1 0 0; 0 0 -1]*v_ned;

	R_enu_ecef = [-sin(lon), cos(lon), 0; 
				-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat); 
				cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]';

	v_ecef = R_enu_ecef*v_enu;
end