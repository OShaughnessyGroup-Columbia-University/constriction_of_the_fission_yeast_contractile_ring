% delmyp7 plot
clear t_del_myp
load('tens_delmyp7_lapv1.mat')
for i = 1:length(tens_cell)
t_del_myp(i, :) = tens_cell{i}(1:22);
end
figure
hold on
t_list = 1:22
r_ring_t = 1.85 - 0.07*(tlist)/2;
phi = acosd(r_ring_t / 1.85);
errorbar(phi(phi< 50), mean(t_del_myp(:, phi< 50)), std(t_del_myp(:, phi< 50)), 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k')
ylim([0, 1000])
tens_exp = [234, 453, 336, 292, 482, 526]; % del myp2
tens_exp_err = [270, 686, 277, 328, 226, 350]/2; % del myp2
phi_exp = 10:10:60; % del myp2
hold on
errorbar(phi_exp,tens_exp,tens_exp_err,'bo', 'Linewidth', 2, 'MarkerFaceColor', 'b');
xlim([0, 70])
ylabel('Ring tension (pN)')
xlabel('Ring position, \phi (deg)')
legend({'Simulation', 'Experiment'}, 'EdgeColor', 'w', 'Location', 'NW')
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
set(gca, 'Position',[0.2 0.2 0.7 0.7]);
box on