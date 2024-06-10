% Function to animate the attitude of a rocket. Based on code developed by
% Alaric Gregoire. 
% 
% Author: Griffin Jourda 8/03/2023
% Inputs
%	t			:	vector of times (s) 
%	q			:	matrix of attitude quaternions (last part being the 
%					scalar component)
%	pauseTime	:	amount of time to pause between each frame update 
%	skipRate	:	amount of data points to skip between each update
function [] = animateAttitude(t, q, pauseTime, skipRate)
	figure(4)
	clf(4)
	ax = gca;
	ax.XLim = [-2, 2];
	ax.YLim = [-2, 2];
	ax.ZLim = [-2, 2];
	ax.Projection = "perspective";
	ax.DataAspectRatio = [1 1 1];
	xlabel('X');
	ylabel('Y');
	zlabel('Z');

	hold on
	grid on

	hgt = hgtransform('Parent',ax);
	hold(ax,'on')

	S = makePrism(0.1, 0.1, 1, q(1, :)', ax);
	view(45, 10)
	set(S,'Parent',hgt)
	
	plot3(ax,[0 1.5],[0 0],[0 0],'r','LineWidth',1.5);
	plot3(ax,[0 0],[0 1.5],[0 0],'g','LineWidth',1.5);
	plot3(ax,[0 0],[0 0],[0 1.5],'b','LineWidth',1.5);

	ts = title(sprintf('Simulation Time %8.2f',t(1)),'FontSize',14);

	for i = 2:skipRate:length(q)
		pause(pauseTime)
		
		qi = q(i, :)';
		ang = 2*acos(qi(4));
		ax = qi(1:3)/norm(qi(1:3));
		R = makehgtform("axisrotate", ax, ang);
		set(hgt, 'Matrix', R);

	    ts.String = sprintf('Simulation Time %8.2f',t(i));
	end

end

function [S] = makePrism(x, y, z, q, ax)
	R = q2R(q);
	xyz = eye(3)*max([x y z]);
	xyz = xyz*R;
	
	c = [0.6 0.6 0.6];
	xs = [0 0 0 0; x x x x;
      	0 x x 0; 0 x x 0;
      	0 x x 0; 0 x x 0];
	ys = [0 0 y y; 0 0 y y;
      	0 0 0 0; y y y y;
      	0 0 y y; 0 0 y y];
	zs = [0 z z 0; 0 z z 0;
      	0 0 z z; 0 0 z z;
      	0 0 0 0; z z z z];
	xs = xs-x/2;
	ys = ys-y/2;
	zs = zs-z/2;
	
	% Body faces
	S(1) = fill3(ax,xs(1,:), ys(1,:), zs(1,:), c);
	S(2) = fill3(ax,xs(2,:), ys(2,:), zs(2,:), c);
	S(3) = fill3(ax,xs(3,:), ys(3,:), zs(3,:), c);
	S(4) = fill3(ax,xs(4,:), ys(4,:), zs(4,:), c);
	S(5) = fill3(ax,xs(5,:), ys(5,:), zs(5,:), c);
	S(6) = fill3(ax,xs(6,:), ys(6,:), zs(6,:), c);
	
	% Body axes
	S(7) = plot3(ax,[0 xyz(1,1)],[0 xyz(2,1)],[0 xyz(3,1)],'r--','LineWidth',1.5);
	S(8) = plot3(ax,[0 xyz(1,2)],[0 xyz(2,2)],[0 xyz(3,2)],'g--','LineWidth',1.5);
	S(9) = plot3(ax,[0 xyz(1,3)],[0 xyz(2,3)],[0 xyz(3,3)],'b--','LineWidth',1.5);
	
	% Fins
	S(10) = fill3(ax, [x/2, 3*x/2, 3*x/2, x/2], [0, 0, 0, 0], [-z/2, -z/2, -z/2 + x, -z/2 + x], 'r');
	S(11) = fill3(ax, -[x/2, 3*x/2, 3*x/2, x/2], [0, 0, 0, 0], [-z/2, -z/2, -z/2 + x, -z/2 + x], 'b');
	S(12) = fill3(ax, [0, 0, 0, 0], [y/2, 3*y/2, 3*y/2, y/2], [-z/2, -z/2, -z/2 + x, -z/2 + x], 'g');
	S(13) = fill3(ax, [0, 0, 0, 0], -[y/2, 3*y/2, 3*y/2, y/2], [-z/2, -z/2, -z/2 + x, -z/2 + x], 'y');
end