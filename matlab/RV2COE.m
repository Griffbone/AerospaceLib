% Convert inertial position and velocity to orbital elements
%
% Author: Griffin Jourda 2/14/2022
%
%	Inputs
% 		rvec	:	ECI position vector
% 		vvec	:	ECI velocity vector
% 		
% 	Outputs
% 		a		:	semimajor axis 
% 		e		:	eccentricity 
% 		i		:	inclination (deg) 
% 		raan	:	right ascension of the ascending node (rad) 
% 		w		:	argument of periapsis (rad)
% 		ta		:	true anomaly (rad) 
% 		arglat	:	argument of latitude (rad) 
% 		truelon	:	true longitude (rad) 
% 		lonper	:	longitude of periapsis (rad)
% 	
% 	Special cases
%		Radial orbit			:	all elements undefined
%		Circular equatorial		:
%		Circular inclined		:
%		Elliptical equatorial	:

function [a, e, i, raan, w, ta, arglat, truelon, lonper] = RV2COE(rvec, vvec, mu)
	tol = 1e-6;
	
	r = norm(rvec); 
	v = norm(vvec); 
	E = v^2/2 - mu/r;
	
	hvec = cross(rvec, vvec); 
	h = norm(hvec);
	
	% Radial orbit 
	if h < tol
		a = nan; 
		e = nan; 
		i = nan; 
		raan = nan; 
		w = nan; 
		ta = nan; 
		arglat = nan; 
		truelon = nan; 
		lonper = nan; 
		return
	end

	nvec = cross([0 0 1], hvec); 
	n = norm(nvec); 
	
	evec = (1/mu)*((v^2 - mu/r)*rvec - dot(rvec, vvec)*vvec);
	e = norm(evec);
	
	i = acos(hvec(3)/h); 
	
	% Semimajor axis 
	if abs(E) > tol
		% elliptical/hyperbolic orbit
		a = -mu/(2*E);
	else
		% parabolic orbit
		a = inf;
	end
	
	% Longitude of ascending node 
	if n > tol
		% inclined orbit case
		raan = acos(nvec(1)/n);
		
		if nvec(2) < 0
			raan = 2*pi - raan;
		end 
	else
		% equatorial orbit case
		raan = nan;
	end
	
	% Argument of periapsis
	if (e > tol) && (i > tol)
		% Eccentric inclined case 
		w = acos(dot(nvec, evec)/(n*e));
		
		if evec(3) < 0
			w = 2*pi - w; 
		end 
	else
		% Circular or non-inclined case 
		w = nan;
	end
	
	% True anomaly 
	if e > tol 
		% Eccentric orbit case
		ta = acos(dot(evec, rvec)/(e*r));
		
		if dot(rvec, vvec) < 0
			ta = 2*pi - ta;
		end
	else
		% Circular orbit case
		ta = nan;
	end
	
	% Argument of laitude (inclined orbit) 
	if n > tol
		arglat = acos(dot(nvec, rvec)/(n*r)); 
		
		if rvec(3) < 0
			arglat = 2*pi - arglat; 
		end
	else 
		arglat = nan;
	end 
	
	% True longitude (circular equatorial) 
	if (e <= tol) && (n <= tol)
		truelon = acos(rvec(1)/r);
		
		if (rvec(2) <  0) || (i > 90)
			truelon = 2*pi - truelon;
		end 
	else
		truelon = nan; 
	end
	
	% True longitude of periapsis (elliptical equatorial) 
	if (e > tol) && (n < tol)
		lonper = acos(evec(1)/e); 
		
		if evec(2) < 0
			lonper = 2*pi - lonper; 
		end 
	else 
		lonper = nan;
	end
end