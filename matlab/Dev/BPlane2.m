function [BR, BT, rp] = BPlane2(vinf1, vinf2, mu)
	S = vinf1/norm(vinf1); 
	h = cross(vinf1, vinf2)/norm(cross(vinf1, vinf2));
	
	B = cross(S, h);
	T = cross(S, [0, 0, 1]); 
	R = cross(S, T);
	
	cosdel = dot(vinf1, vinf2)/(norm(vinf1)*norm(vinf2));
	del = acos(cosdel);
	rp = (mu/norm(vinf1)^2)*((1/cos((pi - del)/2)) - 1);
	
	B = B.*(mu/norm(vinf1)^2)*sqrt((1 + ((norm(vinf1)^2)*rp)/mu)^2 - 1);

	BR = dot(B, R); 
	BT = dot(B, T); 
end