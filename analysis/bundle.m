function bundle = bundle(rbead)
% load('sizep5_8_10min.mat')
% clear rm zm
% rbead(:,ifor) = [];
[phiact, ract, zact] = cart2pol(rbead(1,:), rbead(2,:), rbead(3,:));
rb = .05;

phiact = mod(phiact, 2*pi);
edges = 0:0.0250:max(ract);
n = 50;
% figure
% hold all
slice = zeros(n,length(phiact));
bundle = zeros(n,length(phiact));
% clear slice
for i=1:n
    phi0 = (i-1)*2*pi/n;
    phi1 = i*2*pi/n;
    slice(i,:) = (phiact < phi1 & phiact > phi0).';
    cnts = histcounts(ract(logical(slice(i,:))), edges);
    [val, idx] = max(cnts);
    rm(i) = (edges(idx)+edges(idx+1))/2;
end
rm = movmean(rm, 5);
for i=1:n
    bundle(i,:) = (logical(slice(i,:)) & ...
                 ract < rm(i) + rb & ...
                 rm(i) - rb < ract);
end
sum(sum(bundle))/length(rbead);

% figure
% hold on
% scatter(phiact, ract)
% plot(2*pi/n*(1:n), rm)
% plot(2*pi/n*(1:n), rm + rb)
% plot(2*pi/n*(1:n), rm - rb)
% ylim([0.5, 1.25])

% [phiact, ract, zact] = cart2pol(rbead(1,:), rbead(2,:), rbead(3,:));
% phiact = mod(phiact, 2*pi);
edges = -0.5:0.005:0.5;
slice = zeros(n,length(phiact));
for i=1:n
    phi0 = (i-1)*2*pi/n;
    phi1 = i*2*pi/n;
    slice(i,:) = (phiact < phi1 & phiact > phi0).';
    cnts = histcounts(zact(logical(slice(i,:))), edges);
%     plot(cnts)
    [val, idx] = max(cnts);
    zm(i) = (edges(idx)+edges(idx+1))/2;
end
zm = movmean(zm, 5);

% figure
% hold on
% scatter(phiact, zact)
% plot(2*pi/n*(1:n), movmean(zm, 5))
% plot(2*pi/n*(1:n), movmean(zm, 5)+rb)
% plot(2*pi/n*(1:n), movmean(zm, 5)-rb)
% ylim([-0.25, 0.25])

for i=1:n
    bundle(i,:) = (logical(bundle(i,:)) & logical(slice(i,:)) & ...
                 zact < zm(i) + rb & ...
                 zm(i) - rb < zact);
end
bundle = sum(bundle,1);