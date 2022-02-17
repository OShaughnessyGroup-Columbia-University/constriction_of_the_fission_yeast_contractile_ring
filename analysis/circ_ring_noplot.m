% get the center of all beads
center = mean(rbead,2);

% recenter rbead at the center
r = [];
for i = 1:3
    r = [r; rbead(i,:) - center(i)];
end

% get positions of all beads in cylindrical coordinates
[th, rho, z] = cart2pol(r(1,:),r(2,:),r(3,:));
[~,ind] = sort(th);
th = th(ind);
rho = rho(ind);
z = z(ind);

frho = fit(th', rho', 'smoothingspline','SmoothingParam', 1-1e-2);
fz = fit(th', z', 'smoothingspline','SmoothingParam', 1-1e-2);

[x,y,z] = pol2cart(th', frho(th'), fz(th'));

rbeadstore = rbead;
rmyostore = rmyo;
rbead = r;
rmyo(1,:) = rmyo(1,:) - center(1);
rmyo(2,:) = rmyo(2,:) - center(2);
rmyo(3,:) = rmyo(3,:) - center(3);
rbead = rbeadstore;
rmyo = rmyostore;

%% calculate the length of the ring
% circularly shift x, y and z
xc = circshift(x,1);
yc = circshift(y,1);
zc = circshift(z,1);
% calculate the difference
xd = x - xc;
yd = y - yc;
zd = z - zc;
% length of each segment
lseg = sqrt(xd.^2 + yd.^2 + zd.^2);
% circumference of the ring
circring = sum(lseg);