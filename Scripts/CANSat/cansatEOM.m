% Wrapper of attitude dynamics with the ability to propagate some
% additional things useful for the cansat
function [ydot] = cansatEOM(t, y, I, controlfun)
	q = y(1:4);
	w = y(5:7);
	b = y(8:10);
	b = b/norm(b);

	[ddot] = attitudeEOM(t, y, I, controlfun);
	bdot = -cross(w, b);
	
	ydot = [ddot; bdot];
end