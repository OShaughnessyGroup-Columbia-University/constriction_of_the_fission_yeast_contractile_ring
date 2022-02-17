
% create bins
xbin = linspace(0,circring,100);
tsum = zeros(1,99);
% load file
    tension_noplot
    for j = 1:length(xbin)-1
        if any(and(xfit>xbin(j), xfit<xbin(j+1)))
            tsum(j) = tsum(j) + mean(tens(and(xfit>xbin(j), xfit<xbin(j+1))));
        end
    end
    clear xfit
    clear tens

plot(1:length(tsum),tsum,'o-b')
figure
% get the center of all beads
center = mean(rbead,2);

% recenter rbead at the center
r = [];
for i = 1:3
    r = [r; rbead(i,:) - center(i)];
end

% get positions of all myosin in cylindrical coordinates
[th, rho, z] = cart2pol(r(1,:),r(2,:),r(3,:));
[~,ind] = sort(th);
th = th(ind);
rho = rho(ind);
z = z(ind);

frho = fit(th', rho', 'smoothingspline','SmoothingParam', 1-1e-2);
fz = fit(th', z', 'smoothingspline','SmoothingParam', 1-1e-2);

[x,y,z] = pol2cart(th', frho(th'), fz(th'));

% circularly shift x, y and z
xc = circshift(x,1);
yc = circshift(y,1);
% calculate the difference
xd = x - xc;
yd = y - yc;
% length of each segment
lseg = sqrt(xd.^2 + yd.^2);
% position (horizontal axis of T-x plot)
xfit = cumsum(lseg);

% generate a matrix for the color map
nside = 2000;
cmap = nan(nside);
% generate scale bars for x and y
xscale = linspace(min(y)-.1,max(y)+.1,nside);
yscale = xscale;
% width of the ring
rwidth = .095;

% give each pixel on the color map the value of tension
for i = 1:nside
    for j = 1:nside
        dx = xscale(i)-x;
        dy = yscale(j)-y;
        dist2 = dx.*dx + dy.*dy;
        [mindist2, ind] = min(dist2);
        if mindist2 < rwidth^2
            ind = round(xfit(ind)/circring*100);
            if ind == 0
                ind = 1;
            elseif ind >99
                ind = 99;
            end
            cmap(i,j) = tsum(ind);
        end
    end
end
colormap('hot')
imagescwithnan(cmap,jet,[1 1 1])
axis equal