eclipse = []; 

for i = 1:size(t)
	jdi = jd0 + t(i)/86400;
	
	s = sunVector(jdi);
	r = y(i, 1:3)';

	eclipse = [eclipse, determineEclipse(r, s)];
end

plot(t/3600, eclipse)


function [shadow] = shadow(r, s)
	alf_umb = 0.264121687*pi/180; 
	alf_pen = 0.269007205*pi/180;
	shadow = 0;

	if dot(r, s) < 0
		zeta = acosd(dot(-s, r)/(norm(-s)*norm(r)));

		horiz = norm(r)*cos(zeta);
		vert = norm(r)*sin(zeta); 
		x = 6378000/sin(alf_pen);

		penvert = tan(alf_pen)*(x + horiz);

		if vert < penvert 
			shadow = 1;
			y = 6378000/sin(alf_umb); 
			umbvert = tan(alf_umb)*(y - horiz);

			if vert < umbvert
				shadow = 2;
			end
		end
	end
end 













function [eclipse] = determineEclipse(r, s) 
	theta = asin(6378000/norm(r));
	psi = acos(dot(r, s)/(norm(r)*norm(s)));

	if psi < theta 
		eclipse = true;
	else 
		eclipse = false;
	end
end