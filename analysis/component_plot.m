% component plot
% load('delmyp7_3_2min.mat');
clear rcdc15 r qcdc15 zcdc15
%% actin: grey
figure
axis square
hold all
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
plot(rbead(1,1:ipt(1)), rbead(2,1:ipt(1)), 'Color', [.4 .4 .4])
for i=1:length(ipt)
    plot(rbead(1, 1+ifor(i):ipt(i)), rbead(2, 1+ifor(i):ipt(i)), 'Color', [.4 .4 .4])
end
%% myp2: green
figure
axis square
hold all
if(~isempty(rmyp))
    plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
    scatter(rmyp(1,:),rmyp(2,:),rcapmyp*1000,'Green','filled', 'MarkerFaceAlpha', 0.3)
end
%% myo2: red
% figure
% axis square
% hold all
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
scatter(rmyo(1,:),rmyo(2,:),rcapmyo_short*1000,'Red','filled', 'MarkerFaceAlpha', 0.4)
%% CDC15: brown, 54 nm out from myo2
figure
axis square
hold all
plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi),'black')
[qcdc15, rcdc15, zcdc15] = cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
rcdc15 = rcdc15 + 0.05;
[rcdc15(1,:), rcdc15(2,:), rcdc15(3,:)] = pol2cart(qcdc15, rcdc15, zcdc15);
scatter(rcdc15(1,:), rcdc15(2,:), 50, [.8 .4 0], 'filled')
