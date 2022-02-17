% component plot
% load('delmyp7_3_2min.mat');
clear rcdc15 r qcdc15 zcdc15
%% actin: grey
figure
axis square
hold all
xlim([-1.1*r_ring, 1.1*r_ring])
ylim([-1.1*r_ring, 1.1*r_ring])
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
[~, rract] = cart2pol(rbead(1,:), rbead(2,:));
% rbead = rbead(:, rract < r_ring - 0.05 & rract > r_ring > 0.175);
% iwhisk = find(0.03 > r_ring - rract | 0.18 < r_ring - rract);
iwhisk = find(0.0 > r_ring - rract | 1e6 < r_ring - rract);
plot(rbead(1,1:ipt(1)), rbead(2,1:ipt(1)), 'Color', [.4 .4 .4])
plot([.8*r_ring, .8*r_ring + 0.2], [.8*r_ring, .8*r_ring], 'k', 'LineWidth', 3)
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
% r95-r05;
%% myp2: green
% figure
axis square
hold all
t=0:0.1:2*pi;
xlim([-1.1*r_ring, 1.1*r_ring])
ylim([-1.1*r_ring, 1.1*r_ring])
if(~isempty(rmyp))
    [~, rrmyp] = cart2pol(rmyp(1,:), rmyp(2,:));
%     iwhisk = (0.03 > r_ring - rrmyp | 0.17 < r_ring - rrmyp);
    iwhisk = (0.0 > r_ring - rrmyp | 10^6 < r_ring - rrmyp);
    plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
    plot([.8*r_ring, .8*r_ring + 0.2], [.8*r_ring, .8*r_ring], 'k', 'LineWidth', 3)
%     scatter(rmyp(1,:),rmyp(2,:),rcapmyp*1003,'Green','filled', 'MarkerFaceAlpha', 0.3)
    for i=1:length(rmyp)
        if(~iwhisk(i))
            xc = rmyp(1,i) + rcapmyp*sin(t);
            yc = rmyp(2,i) + rcapmyp*cos(t);
            fill(xc, yc, [0, 0.6, 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
        end
    end
end
%% myo2: red
% figure
axis square
hold all
xlim([-1.1*r_ring, 1.1*r_ring])
ylim([-1.1*r_ring, 1.1*r_ring])
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
plot([.8*r_ring, .8*r_ring + 0.2], [.8*r_ring, .8*r_ring], 'k', 'LineWidth', 3)
% scatter(rmyo(1,:),rmyo(2,:),rcapmyo_short*1000,'Red','filled', 'MarkerFaceAlpha', 0.4)
t=0:0.1:2*pi;
for i=1:length(rmyo)
    [q, ~] = cart2pol(rmyo(1,i), rmyo(2,i));
    qhat = [-sin(q), cos(q)];
    rhat = [cos(q), sin(q)];
    xc = rmyo(1,i) + rcapmyo_short*sin(t)*cos(q) - rcapmyo_long*cos(t)*sin(q);
    yc = rmyo(2,i) + rcapmyo_long*cos(t)*cos(q) + rcapmyo_short*sin(t)*sin(q);
    fill(xc, yc, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
end
%% CDC15: brown, 54 nm out from myo2
% figure
axis square
hold all
xlim([-1.1*r_ring, 1.1*r_ring])
ylim([-1.1*r_ring, 1.1*r_ring])
t=0:0.1:2*pi;
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
[qcdc15, rcdc15, zcdc15] = cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
rcdc15 = rcdc15 + 0.05;
[rcdc15(1,:), rcdc15(2,:), rcdc15(3,:)] = pol2cart(qcdc15, rcdc15, zcdc15);
% scatter(rcdc15(1,:), rcdc15(2,:), 50, [.8 .4 0], 'filled')
for i=1:length(rcdc15)
    [q, ~] = cart2pol(rcdc15(1,i), rcdc15(2,i));
    qhat = [-sin(q), cos(q)];
    rhat = [cos(q), sin(q)];
    xc = rcdc15(1,i) + .035*sin(t);
    yc = rcdc15(2,i) + .035*cos(t);
    fill(xc, yc, [.85, .65, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
end
plot([.8*r_ring, .8*r_ring + 0.2], [.8*r_ring, .8*r_ring], 'k', 'LineWidth', 3)