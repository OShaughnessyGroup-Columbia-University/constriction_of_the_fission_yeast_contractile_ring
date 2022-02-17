% plots cross-section of ring and performs some calculations.
figure
% plot myp c.s. and figure out where the ring lies.
if exist('rmyp','var') & ~isempty(rmyp)
	[theta,r,z]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));

	for k=1:length(r)
		xc = z(k);
		yc = r_ring-r(k);
		rb = 0.1;
		x = rb*sin(-pi:0.1*pi:pi) + xc;
		y = rb*cos(-pi:0.1*pi:pi) + yc;
		c = [32 180 32]./255;
		fill(x, y, c, 'EdgeColor','None')
		hold on;
	end
	
	rmypMedian = r_ring - median(r);
	zmypMedian = median(z);
		
	xc = zmypMedian;
	yc = rmypMedian;
	rmyx = 0.1;
	rb = 0.1;
	x = rb*sin(-pi:0.01*pi:pi) + xc;
	y = rb*cos(-pi:0.01*pi:pi) + yc;
	c = [1 1 1];
	fill(x, y, c,'FaceColor','None','EdgeColor','k', 'LineWidth', 4,...
        'LineStyle', '--')
end

% plot myo c.s.
if exist('rmyo','var') & ~isempty(rmyo)
	[theta,r,z]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));

	for k=1:length(r)
		xc = z(k);
		yc = r_ring-r(k);
		rb = 0.051;
		x = rb*sin(-pi:0.1*pi:pi) + xc;
		y = rb*cos(-pi:0.1*pi:pi) + yc;
% 		c = [1 0 0];
		c = [164 25 25]./255;
		fill(x, y, c, 'EdgeColor','None')
		hold on;
	end
	
	rmyxMedian = r_ring - median(r);
	zmyxMedian = median(z);
	
	xc = zmyxMedian;
	yc = rmyxMedian;
	rmyx = 0.051;
	rb = 0.051;
	x = rb*sin(-pi:0.01*pi:pi) + xc;
	y = rb*cos(-pi:0.01*pi:pi) + yc;
	ck = [1 1 1];
	
	if ~exist('rmyp','var')
		rmyp = [];
    end
	
    xmyo = x;
    ymyo = y;
	
end

% plot actin c.s. and identify bundle as a rectangle.
if exist('rbead','var')
	flgB = [];
	[thetaB,rB,zB]=cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));

	for k=1:length(rB)
		xc = zB(k);
		yc = r_ring-rB(k);
		rb = 0.0025;
		
		%d1sq = (zmypMedian - xc)^2 + (rmypMedian - yc)^2;
		
		x = rb*sin(-pi:0.1*pi:pi) + xc;
		y = rb*cos(-pi:0.1*pi:pi) + yc;
		%c = [0 0 1];
		c = [0.5 0.5 0.5];
		fill(x, y, c, 'EdgeColor','None')
		hold on;
	end
	
	c = [0 0.5 0];
end

if isempty(rmyp)
    fill(xmyo, ymyo, ck,'FaceColor','None','EdgeColor','k','LineWidth', 4)
end
    
axis equal;
xlim([-0.15 0.15]);
ylim([0 0.4]);
%xlim([-1 1]);
%ylim([0 2]);
%xlim([-0.4 0.4]);
%ylim([0 0.7]);
view(180,90); 
axis off;