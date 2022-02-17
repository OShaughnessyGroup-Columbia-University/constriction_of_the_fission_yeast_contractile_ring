% component plot
% load('delmyp7_3_2min.mat');
clear rcdc15 r qcdc15 zcdc15
%% for all
figure
axis square
hold all
xlim([-1.1*r_ring, 1.1*r_ring])
ylim([-1.1*r_ring, 1.1*r_ring])
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
%% actin: grey
[~, rract] = cart2pol(rbead(1,:), rbead(2,:));
% rbead = rbead(:, rract < r_ring - 0.05 & rract > r_ring > 0.175);
iwhisk = find(0.0 > r_ring - rract & 0 < r_ring - rract);
plot(rbead(1,1:ipt(1)), rbead(2,1:ipt(1)), 'Color', [.4 .4 .4])
for i=1:length(ipt)
    filend = find(1+ifor(i) < iwhisk & iwhisk < ipt(i), 1);
    if(isempty(filend))
        filend = ipt(i);
    else
        filend = iwhisk(filend);
    end
    plot(rbead(1, 1+ifor(i):filend), rbead(2, 1+ifor(i):filend), 'Color', [.4 .4 .4])
%     plot(rbead(1, 1+ifor(i):ipt(i)), rbead(2, 1+ifor(i):ipt(i)), 'Color', [.4 .4 .4])
end
r95 = quantile(rract(~(.05 > r_ring - rract | 0.18 < r_ring - rract)),0.95);
r05 = quantile(rract(~(.05 > r_ring - rract | 0.18 < r_ring - rract)),0.05);
rwhisk = length(iwhisk)/length(rract)
plot([.8*r_ring, .8*r_ring + 0.1], [.8*r_ring, .8*r_ring], 'k', 'LineWidth', 3)
r95-r05;