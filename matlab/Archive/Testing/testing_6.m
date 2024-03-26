clc; clear; 
addpath('C:\Users\Griffin\OneDrive\GaTech\AE 4361\Library')

% Spacraft properties
lx = 3*10/100; 
ly = 2*10/100; 
lz = 1*10/100; 
m = 7; 

I = [(1/12)*m*(ly^2 + lz^2), 0, 0;
	0, (1/12)*m*(lx^2 + lz^2), 0;
	0, 0, (1/12)*m*(lx^2 + ly^2)];	% kg*m^2

% Initial conditions
rp = 600e3 + 6371e3; % 6378e3; 
ra = 39700e3 + 6371e3; % 10*6378e3; 
a = (rp + ra)/2; 
e = ra/a - 1;
T = 2*pi*sqrt(a^3/398600435507000);

[r0, v0] = COE2RV(a, e, deg2rad(63.4), deg2rad(240), deg2rad(270), deg2rad(180), 398600435507000);
w0 = deg2rad([0; 0; 0]);
q0 = [0; 0; 0; 1];
y0 = [r0'; v0'; w0; q0];

tspan = 0:0.1:T*2;
y = ode4(@(t, y) fullStateEOM(t, y, I), tspan, y0); 

%% Plotting 
skip = 500; 
tds = tspan(1:skip:end); 
yds = y(1:skip:end, :);
animateOrbit(tds, yds, 100000);
% animateAttitude(tspan, y(:, 7:end), 100, false, I)

%% Functions
function [ydot] = mTorque(t, y, I, kp, kd)
	w = y(1:3); 
	q = y(4:7); 

	% Magnetic field 
% 	B = [30000e-9; 0; 0];	% inertial frame 
% 	B = q2R(q)*B;			% body frame 
% 	b = B/norm(B);			% unit field vector 
	
% 	m = (0.075)*cross(w/norm(w), b);		% applied dipole 
% 	M = cross(m, B*10);						% magnetorquer torque

	% Control 
	if t < 15
		qc = [0; 0; 0; 1];
	elseif t > 15
		qc = v2q(deg2rad([90; 90; 0]));
	end

	dq = qMult(q, q_inv(qc));
	M = -kp*dq(1:3) - kd*w;

	% Attitude dynamics
	wdot = I\(M - cross(w, I*w));
	qdot = (1/2)*qMult([w; 0], q); 
	ydot = [wdot; qdot];
end


function [ydot] = fullStateEOM(t, y, I)
	r = y(1:3); 
	v = y(4:6); 
	w = y(7:9); 
	q = y(10:13); 

	% Orbital Dynamics 
	a = -(398600435507000/norm(r)^3)*r;

	% Attitude control (nadir-pointing) 
% 	kp = 0.1; kd = 0.1; 
% 	xhat = -r/norm(r); 
% 	yhat = cross(r, v)/norm(cross(r, v));
% 	zhat = cross(xhat, yhat);
% 	R = [xhat, yhat, zhat]'; 
% 
% 	phi = R2v(R);
% 	qc = v2q(phi);
% 	dq = qMult(q, q_inv(qc));
% 	M = -kp*dq(1:3) - kd*w;
	M = [0; 0; 0];

	% Attitude dynamics 
	wdot = I\(M - cross(w, I*w));
	qdot = (1/2)*qMult([w; 0], q); 

	% Create state vector derivative
	ydot = [v; a; wdot; qdot];
end