clear rm zm
rb = rbead;
rb(:,ifor) = [];
[phiact, ract, zact] = cart2pol(rbead(1,:), rbead(2,:), rbead(3,:));
[phimyo, rrmyo, zmyo] = cart2pol(rmyo(1,:), rmyo(2,:), rmyo(3,:));
rb = .07;

phiact = mod(phiact, 2*pi);
phimyo = mod(phimyo, 2*pi);
edges = 0:0.020:max(ract);
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
%     if i>1
%         [csrt, sidx] = sort(cnts, 'descend');
%         [mmc, idx] = min(abs(csrt(1:5) - rm(i-1)));
%         sidx = sidx(idx);
%         length(sidx)
%     else
%         [val, sidx] = max(cnts);
%     end
    [val, sidx] = max(cnts);
    rm(i) = (edges(sidx)+edges(sidx+1))/2;
end
rm = movmean(rm, 5);
for i=1:n
    bundle(i,:) = (logical(slice(i,:)) & ...
                 ract < rm(i) + rb & ...
                 rm(i) - rb < ract);
end
sum(sum(bundle))/length(rbead);

figure
hold on
scatter(phiact, ract)
scatter(phimyo, rrmyo)
plot(2*pi/n*(1:n), rm)
plot(2*pi/n*(1:n), rm + rb)
plot(2*pi/n*(1:n), rm - rb)
ylim([r_ring/2, r_ring+0.1])
mean(rm+rb);
mean(rm-rb);

% [phiact, ract, zact] = cart2pol(rbead(1,:), rbead(2,:), rbead(3,:));
% phiact = mod(phiact, 2*pi);
edges = -0.5:0.05:0.5;
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

figure
hold on
scatter(phiact, zact)
scatter(phimyo, zmyo)
plot(2*pi/n*(1:n), movmean(zm, 5))
plot(2*pi/n*(1:n), movmean(zm, 5)+rb)
plot(2*pi/n*(1:n), movmean(zm, 5)-rb)
ylim([-0.25, 0.25])

for i=1:n
    bundle(i,:) = (logical(bundle(i,:)) & logical(slice(i,:)) & ...
                 zact < zm(i) + rb & ...
                 zm(i) - rb < zact);
end
bundle = sum(bundle,1);
rwhisk = 1-length(find(bundle))/length(bundle)