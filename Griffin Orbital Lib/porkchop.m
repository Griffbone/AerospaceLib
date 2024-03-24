% Function to create data for a porkchop plot of transfers between two
% planets 
% 
% Author: Griffin Jourda 10/12/22
% 
%	Inputs 
%		dp		:	departure planet (1 = Mercury, 2 = Venus, ... 8 = Neptune)
%		ap		:	arrival planet, same as above 
%		jd1d	:	first departure Julian date to consider 
%		jd2d	:	last departure Julian date to consider 
%		jd1a	:	first arrival Julian date to consider 
%		jd2a	:	last arrival Julian date to consider
%		tmin	:	minimum transfer time (days) 
%		tmax	:	maximum transfer time (days)
%		res		:	resolution of the porkchop plot
% 
%	Outputs 
%		vinfd	:	array of departure excess velocities (km/s)
%		vinfa	:	array of arrival excess velocities (km/s)
%		xx
%		yy

function [vinfd, vinfa, xx, yy] = porkchop(dp, ap, jd1d, jd2d, jd1a, jd2a, res)
	load('constants.mat', 'mu_s');
	
	jdd = linspace(jd1d, jd2d, res); 
	jda = linspace(jd1a, jd2a, res); 
	
	vinfd = zeros(res);
	vinfa = zeros(res); 
	
	for i = 1:res % res:-1:1
		jdai = jda(i);
		for j = 1:res
			jddj = jdd(j); 
			
			[rd, vd] = planetState(jddj, dp); 
			[ra, va] = planetState(jdai, ap); 
			[v1, v2] = lambertUV(rd, ra, (jdai - jddj)*86400, mu_s, 1); 
			
			vinfd(i, j) = norm(v1 - vd); 
			vinfa(i, j) = norm(v2 - va);			
		end
	end
	
	[xx, yy] = meshgrid(jdd, jda);
end