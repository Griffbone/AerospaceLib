clc; clear; 
addpath('../../matlab/')

%% Determine MMOIs
m = 350/1000; 
l = 115/1000;
r = 33/1000;

Ixx = 0.5*m*r^2; 
Iyy = 0.25*m*r^2 + (1/12)*m*l^2;
I = [Ixx, 0, 0;
     0, Iyy, 0; 
     0, 0, Iyy];

%% Simulate
q0 = [0; 0; 0; 1];
w0 = [10; 50; -30]*pi/180;
b0 = [1; 0; 0;]; 
y0 = [q0; w0; b0];

tspan = linspace(0, 60*5, 10000);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, y] = ode45(@(t, y) cansatEOM(t, y, I, @(t, y) cansatControl(t, y)), tspan, y0, options);

% Back out command torques
M = zeros(length(tspan), 3);
bvec = zeros(length(tspan), 3);
for i = 1:length(tspan)
	M(i, :) = cansatControl(tspan(i), y(i, :)');
	bvec(i, :) = (q2R(y(i, 1:4)')*[1; 0; 0])';
end

q = y(:, 1:4);
w = y(:, 5:7);
save('sim_data', 'tspan', 'q', 'w');

%% Plot
% animateAttitude(t, y(:, 1:4), 0.01, 10)

figure
subplot(2, 1, 1)
plot(tspan, y(:, 5:7)*180/pi)
xlabel('Time [s]')
ylabel('Body Angular Velocity [deg/s]')

subplot(2, 1, 2)
plot(tspan, M*1000)
xlabel('Time [s]')
ylabel('Torque Commands [mN*m]')

figure 
plot(tspan, y(:, 8:10) - bvec)