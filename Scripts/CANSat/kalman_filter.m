load('sim_data.mat');

Q = eye(9)*1e-4;	% gotta find this for noisy gyro model
R = ones(6, 1)*1e-3;
bias_true = [0; 0; 0]*180/pi;
w_noise = 1e-3;


for i = 1:length(tspan)
	% Predict 
	w_meas =  w(i, :)' + bias_true + randn(3, 1)*w_noise;
	Phi = eye(9) + F(w_meas - bias_post)*dt; 
	x_prior = Phi*x_post;
	P_prior = Phi*P_post*Phi' + Q;
	
	% Update
	H = observation_matrix();
	K = P_prior*H'*(H*P_prior*H' + R);
end

%% Kalman filter functions 
function [H] = observation_matrix()
	I3 = eye(3);
	z3 = zeros(3);
	H = [I3, z3, z3; 
		I3, z3, z3];
end

function [F] = dynamics_matrix(w)
	Omega = -skewSym(w);
	z3 = zeros(3);

	F = [Omega, z3, z3;
		z3, Omega, z3;
		z3, z3, z3];
end