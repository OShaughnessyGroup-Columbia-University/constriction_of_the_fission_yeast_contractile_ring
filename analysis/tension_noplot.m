% % % get the center of all myosin
% % center = mean(rmyo,2);
% % 
% % % r is rbead recentered
% % r = [];
% % for i = 1:3
% %     r = [r; rmyo(i,:) - center(i)];
% % end

r = rmyo;

% get sorted positions of all beads in cylindrical coordinates
[th, rho, z] = cart2pol(r(1,:),r(2,:),r(3,:));
[~,ind] = sort(th);
th = th(ind);
rho = rho(ind);
z = z(ind);

% spline fitting
frho = fit(th', rho', 'smoothingspline','SmoothingParam', 1-1e-2);
fz = fit(th', z', 'smoothingspline','SmoothingParam', 1-1e-2);
fitth = linspace(-pi, pi, 1001);
fitth(end) = [];
fitrho = frho(fitth');
fitz = fz(fitth');

% cartesian coordinates of fitting
[x,y,z] = pol2cart(fitth', fitrho, fitz);

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
% position (horizontal axis of T-x plot)
xfit = cumsum(lseg);
%% Calculate tension
% calculate the tension in each segment of actin filament (between beads)
get_segtens

% shift rbead and rmyo to center
% % % rbc = rbead - repmat(center,1,length(rbead));
rbc = rbead;

% calculate ring tension accross angle fitth
tens = get_tens_ring(segtens,fitth,rbc,ifor,ipt);