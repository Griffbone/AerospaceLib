% Function to compute rotational equations of motion for a rigid body. Body
% angular rate dynamics are governed by Euler's rigid body equations while
% the attitude representation kinematics are governed by quaternion
% kinematics. 
% 
% Author: Griffin Jourda 6/10/24
% 
% Inputs
%		t			:	time (s) 
%		y			:	state vector [q; w] (n/a, rad/s)
%		I			:	moment of inertia tensor (kg*m^2)
%		controlfun	:	attitude controller function - must take time and state 
%						as inputs and return a column vector of moments in N*m
%						about each body axis: M = controlfun(t, y)
%
%	Outputs 
%		ydot		:	state vector derivative [qdot; wdot] (n/a, rad/s^2)

function [ydot] = attitudeEOM(t, y, I, controlfun)
	q = y(1:4);
	q = q/norm(q);	% re-normalize quaternion (floating point precision issues)
	w = y(5:7);

	M = controlfun(t, y);
	wdot = I\(M - cross(w, I*w));
	qdot = (1/2)*qMult([w; 0], q); 

	ydot = [qdot; wdot];
end