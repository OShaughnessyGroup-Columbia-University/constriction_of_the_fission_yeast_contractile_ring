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
		fill(x, y, c, 'FaceAlpha', 0.04, 'EdgeColor','None')
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
	%fill(x, y, c,'FaceColor','None','EdgeColor','k')
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
		c = [164 25 25]./255;
		fill(x, y, c, 'FaceAlpha', 0.04, 'EdgeColor','None')
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
	c = [1 1 1];
	
	if ~exist('rmyp','var')
		rmyp = [];
	end
	
% 	if isempty(rmyp)
% 		fill(x, y, c,'FaceColor','None','EdgeColor','k')
% 	end
	
end

% plot actin c.s. and identify bundle as a rectangle.
if exist('rbead','var')
	flgB = [];
	[thetaB,rB,zB]=cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));

	for k=1:length(rB)
        if(rand < 1)
            xc = zB(k);
            yc = r_ring-rB(k);
            rb = 0.0025;
            r2b = (zmypMedian - xc)^2 + (rmypMedian - yc)^2;

            if r2b <= rmyx^2% & yc >= 0.046
                flgB = [flgB true];
            else
                flgB = [flgB false];
            end

            x = rb*sin(-pi:0.1*pi:pi) + xc;
            y = rb*cos(-pi:0.1*pi:pi) + yc;
            %c = [0 0 1];
            c = [0.5 0.5 0.5];
            fill(x, y, c, 'FaceAlpha', 0.4, 'EdgeColor','None')
%             fill(x, y, c, 'FaceAlpha', 0.3 + min([25*r2b, 0.7]), 'EdgeColor','None')
%             if(abs(xc) < 0.06 && 0.04 < yc && yc < 0.17)
%                 fill(x, y, c, 'FaceAlpha', 0.3, 'EdgeColor','None')
%             else
%                 fill(x, y, c, 'FaceAlpha', 1, 'EdgeColor','None')
%             end
            hold on;
        end
	end
	
	%plot(zB,r_ring-rB,'b.');
	
	wdthBTemp = zB(logical(flgB));
	thckBTemp = rB(logical(flgB));
	
	wdthB = prctile(wdthBTemp,99) - prctile(wdthBTemp,1)
	thckB = prctile(thckBTemp,99) - prctile(thckBTemp,1)
	
	c = [0 0.5 0];
	
	
	left = zmypMedian - wdthB/2;
	right = left + wdthB;
	bottom = rmypMedian - thckB/2;
	top = bottom + thckB;
	x = [left left right right];
	y = [bottom top top bottom];
	
% 	fill(x, y, c,'FaceColor','None','EdgeColor','k')
	
	
	
end

axis equal;
xlim([-0.25 0.25]);
ylim([0 0.4]);
% ylim([0 r_ring]);
% ylim([0 0.4]);
plot(xlim, [0,0], 'k', 'LineWidth', 1)
%xlim([-1 1]);
%ylim([0 2]);
%xlim([-0.4 0.4]);
%ylim([0 0.7]);
view(180,90); 
axis off;

% print percentage of non-whisker subunits.
sum(flgB)/size(rbead,2)*100