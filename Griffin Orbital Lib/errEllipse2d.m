function [] = errEllipse2d(r, cov)
	[u, lam] = eig(cov);
	
	[~, order] = sort(max(lam), 'descend');
	u = u(:, order); 
	lam = lam(:, order);
	
	radii = max(lam);
	a = sqrt(radii(1)); 
	b = sqrt(radii(2)); 
	
	% Generate ellipse
	x = linspace(-a, a, 1000); 
	y = (b/a)*sqrt(a^2 - x.^2);
	x = [x, x(end:-1:1)];
	y = [y, -y];
	
	xp = zeros(1, length(x)); 
	yp = zeros(1, length(x));
	
	for i = 1:length(x)
		p = [x(i); y(i)];
		pp = u*p;
		xp(i) = pp(1);
		yp(i) = pp(2);
	end
		
	% Plot error ellipses
	hold on; grid on; axis equal;
	plot([r(1), u(1, 1) + r(1)], [r(2), u(2, 1) + r(2)], 'r', 'HandleVisibility', 'off')
	plot([r(1), u(1, 2) + r(1)], [r(2), u(2, 2) + r(2)], 'b', 'HandleVisibility', 'off')
	
	plot(xp + r(1), yp + r(2), 'k')
	plot(xp*2 + r(1), yp*2 + r(2), 'k')
	plot(xp*3 + r(1), yp*3 + r(2), 'k')
	
	lpos1 = u*[a; 0] + [r(1); r(2)];
	lpos2 = u*[2*a; 0] + [r(1); r(2)];
	lpos3 = u*[3*a; 0] + [r(1); r(2)];
	text(lpos1(1)*0.95, lpos1(2)*0.95, '1-\sigma')
	text(lpos2(1)*0.95, lpos2(2)*0.95, '2-\sigma')
	text(lpos3(1)*0.95, lpos3(2)*0.95, '3-\sigma')
end