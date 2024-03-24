% Function to find VNB basis vectors expressed in inertial frame 
% 
% Author: Griffin Jourda 10/13/22
%
%	Inputs 
%		r	:	inertial position 
%		v	:	inertial velocity 
% 
%	Outputs 
%		vhat	:	v basis vector 
%		nhat	:	normal basis vector 
%		bhat	:	binormal basis vector 

function [vhat, nhat, bhat] = VNB(r, v) 
	vhat = v/norm(v); 
	nhat = cross(r, v)/norm(cross(r, v)); 
	bhat = cross(vhat, nhat); 
end