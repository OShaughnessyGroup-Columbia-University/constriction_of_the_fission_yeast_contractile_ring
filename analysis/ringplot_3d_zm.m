figure
hold on

for i = 1:length(ifor)
    plot3(rbead(1, ifor(i):ipt(i)), rbead(2, ifor(i):ipt(i)), rbead(3, ifor(i):ipt(i)), 'Color', [0.5, 0.5 0.5])
%     view([10, 10, 10])
end
plot3(rmyo(1, :), rmyo(2, :), rmyo(3, :), 'o', 'Color', 'r', 'MarkerFaceColor', 'r')
plot3(rmyp(1, :), rmyp(2, :), rmyp(3, :), 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'Markersize', 8)
plot3(rbead(1, ifor), rbead(2, ifor), rbead(3, ifor), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'Markersize', 3)
grid on
axis equal
xlim([-2, 2])
ylim([-2, 2])
zlim([-2, 2])

%% show "septum"
% phi = 0:0.02*pi:2*pi;
% z = linspace(-max(abs(rbead(3, :))), max(abs(rbead(3, :))), 100);
% [PHI, Z] = meshgrid(phi, z);
% surf(r_ring*cos(PHI), r_ring*sin(PHI), Z, 'Facealpha', 0.4, 'edgecolor', 'none')