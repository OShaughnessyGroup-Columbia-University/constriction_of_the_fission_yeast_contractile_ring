% plots cross-section of ring and performs some calculations.
figure
% plot myp c.s. and figure out where the ring lies.
% q0 = 0;
dq = 0.100/r_ring;
r1 = 0.020*sqrt(log(X1(irun)/3))/2;
if r1 < 0.001
    r1 = 0.001;
end


if exist('rmyp','var') & ~isempty(rmyp)
    [theta,r,z]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));
    this_x = q0 < theta & theta < q0+dq;
    theta = theta(this_x);
    r = r(this_x);
    z = z(this_x);
    
    for k=1:length(r)
        xc = z(k);
        yc = r_ring-r(k);
        rb = 0.1;
        x = rb*sin(-pi:0.1*pi:pi) + xc;
        y = rb*cos(-pi:0.1*pi:pi) + yc;
        c = [32 180 32]./255;
        viscircles([xc; yc]', 0.1*ones(1, length(xc)), 'Color', c);
        viscircles([xc; yc]', 0.1*ones(1, length(xc))/2, 'Color', c, 'LineStyle', '--');

%         fill(x, y, c, 'FaceAlpha', 0.14, 'EdgeColor','None')
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
end

% plot myo c.s.
if exist('rmyo','var') & ~isempty(rmyo)
    [theta,r,z]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
    this_x = q0 < theta & theta < q0+dq;
    theta = theta(this_x);
    r = r(this_x);
    z = z(this_x);
    
    for k=1:length(r)
        xc = z(k);
        yc = r_ring-r(k);
        rb = 0.051;
        x = rb*sin(-pi:0.1*pi:pi) + xc;
        y = rb*cos(-pi:0.1*pi:pi) + yc;
        % 		c = [1 0 0];
        c = [164 25 25]./255;
        viscircles([xc; yc]', rcapmyo_short*ones(1, length(xc)), 'Color', c);
        viscircles([xc; yc]', rcapmyo_short*ones(1, length(xc))/2, 'Color', c, 'LineStyle', '--');
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
    this_x = q0 < thetaB & thetaB < q0+dq;
    thetaB = thetaB(this_x);
    rB = rB(this_x);
    zB = zB(this_x);
    nfil = length(rB);
    
    for k=1:length(rB)
        if(~any(ifor == k))
            xc = zB(k);
            yc = r_ring-rB(k);
            rb = 0.0025;
            
            %d1sq = (zmypMedian - xc)^2 + (rmypMedian - yc)^2;
            
            x = rb*sin(-pi:0.1*pi:pi) + xc;
            y = rb*cos(-pi:0.1*pi:pi) + yc;
            %c = [0 0 1];
            c = [0.5 0.5 0.5];
            %         fill(x, y, c, 'FaceAlpha', 1, 'EdgeColor','None')
            viscircles([xc; yc]', r1*ones(1, length(xc)), 'Color', c);
            hold on;
        end
    end
    
    c = [0 0.5 0];
end

if isempty(rmyp)
    fill(xmyo, ymyo, ck,'FaceColor','None','EdgeColor','k','LineWidth', 4)
end

axis equal;
xlim([-0.15 0.15]);
ylim([0 0.3]);
%xlim([-1 1]);
%ylim([0 2]);
%xlim([-0.4 0.4]);
%ylim([0 0.7]);
view(180,90);
axis off;