clc; clear; 

mu = 324858.592;													% km^3/s^2
vinf1 = [-3.32760505287096, 3.13387574455063, -0.878312716887848];	% km/s 
vinf2 = [-4.59935748832974, -2.58249103377671, -1.90250836681768];	% km/s

[BR, BT, rp] = BPlane2(vinf1, vinf2, mu);
D = sqrt(BR^2 + BT^2);

eps1 = norm(vinf1)^2/2;
h1 = D*norm(vinf1);
e1 = sqrt(1 + (2*eps1*h1^2)/mu^2);
a1 = -mu/(2*eps1);

eps2 = norm(vinf2)^2/2;
% h2 = 

%% 
clc; clear;
rpplanet = 6051.8;
mu = 324858.592;

D = 20989.022861188;
rp = 11113.8350235628;
vinf1 = 4.65462842757851;
vinf2 = 5.60739577083723;

v0 = [0; vinf1; 0];
r0 = [D; -0.616e6; 0];

% Propagate to periapsis 
t1 = linspace(0, 4*86400, 5000);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 60, 'InitialStep', 60, 'Events', @(t, y) PeEvent(t, y, mu));
[t1, y1, ~, ~] = ode45(@(t, y) twoBodyEOM(t, y, mu), t1, [r0; v0], options); 

% Propagate other asymptote 
t2 = linspace(-max(t1), max(t1), 5000); 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 60, 'InitialStep', 60);
[t2, y2] = ode45(@(t, y) twoBodyEOM(t, y, mu), t2, [y1(end,1:3)'; y1(end,4:6)'], options); 

% Propagate other asymptote (with maneuver) 
v0 = y1(end,4:6) + (y1(end,4:6)/norm(y1(end,4:6)))*1;
[t3, y3] = ode45(@(t, y) twoBodyEOM(t, y, mu), t2, [y1(end,1:3)'; v0'], options);


close all

y = [y1; y2];
y2 = [y1; y3];
m = (y(end,2) - y(end - 50, 2))/(y(end, 1) - y(end - 50, 1)); 
disp(90 + atand(m))
m2 = (y2(end,2) - y2(end - 50, 2))/(y2(end, 1) - y2(end - 50, 1)); 
disp(90 + atand(m2))

plot(y(:, 1)/rpplanet, y(:, 2)/rpplanet) 
axis equal
hold on 

plot(y2(:,1)/rpplanet, y2(:,2)/rpplanet)

x = linspace(-10, 10, 1000);
% plotLine(y(end, 1:2)/rpplanet, m, x)
plotLine(y2(end, 1:2)/rpplanet, m2, x)
plotCircle(1, 0, 0)
% xline(D/rpplanet)
xlim([-10, 10])
ylim([-10, 10]) 
% axis equal

%% Functions 
function [] = plotCircle(r, cx, cy)
	ts = linspace(0, 2*pi, 1000);
	x = r*cos(ts);
	y = r*sin(ts);

	fill(x + cx, y + cy, 'y', 'FaceAlpha', 0.25)
end


function [] = plotLine(p, m, x)
	y = m*(x - p(1)) + p(2);
	plot(x, y, 'k')
end

function [tadot, isterminal, direction] = PeEvent(~, y, mu)
	r = y(1:3); 
	v = y(4:6); 
	
	e = (1/mu)*((norm(v)^2 - mu/norm(r))*r - dot(r, v)*v);
	tadot = abs(1 - dot(e, r)/(norm(e)*norm(r)));
	
	if tadot < 1e-6
		tadot = 0; 
% 		disp(acosd(tadot))
	end 
	
	isterminal = 1; 
	direction = 0;
end