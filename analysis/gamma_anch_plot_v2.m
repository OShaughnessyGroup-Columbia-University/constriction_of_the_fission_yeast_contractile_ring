% figure
% load('gm_2/gm_2_10_600sec.mat')

gm_theory = 0:100:3500;
l = mean(diff(ifor)-1)*0.1;
rho_myo = length(rmyo)/(2*pi*r_ring);
rho_myp = length(rmyp)/(2*pi*r_ring);
gamma_myo = fhead/vmyo0;
gamma_myp = fheadmyp/vmyo0;
vnode_theory = l*(rho_myo*fhead + rho_myp*fheadmyp)*(1 - vpol/vmyo0)./...
    (gm_theory*rho_myo*l/16 + gmsol + l*(2*rho_myo*gamma_myo + rho_myp*gamma_myp));
% vnode_theory = l*(rho_myo*fhead + rho_myp*fheadmyp)*(1 - vpol/vmyo0)./...
%     (gm_theory + gmsol + l*(2*rho_myo*gamma_myo + rho_myp*gamma_myp));
t_theory = 16*l/2*(rho_myo*fhead + rho_myp*fheadmyp)*...
    (1 - (vpol+vnode_theory)/vmyo0);
% tens_analytical = 16*1.3 * (1-vpol_exp/.24) .* (n_myo_exp+n_myp_exp)...
%     ./r_ring_array/2/pi .* l_exp/2;

figure
set(gca, 'Position',[0.2 0.2 0.6 0.7]);
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
box on
hold on
% ylim([0, 0.2])
errorbar(gamma_anch, vnode_mean, vnode_std, 'ko', 'LineWidth', 2, ...
    'MarkerFaceColor', 'auto', 'MarkerSize', 10, 'CapSize', 12)
% plot(gm_theory, vnode_theory, 'LineWidth', 3)
plot(xlim, [.028, .028], 'r--', 'LineWidth', 2)
plot([2e3,2e3], ylim, 'k--', 'LineWidth', 2)
ylabel('vnode (mm/s)')
xlabel('ganch (pN s/mm)')
xlim([0, 3500])

figure
set(gca, 'Position',[0.2 0.2 0.6 0.7]);
set(gca, 'Linewidth', 3)
set(gca, 'FontSize', 30)
box on
hold on
errorbar(gamma_anch, tens_ave, tens_std, 'ko', 'LineWidth', 2, ...
    'MarkerFaceColor', 'auto', 'MarkerSize', 10, 'CapSize', 12)
% plot(gm_theory, t_theory, 'LineWidth', 3)
plot(xlim, [640 640], 'k--', 'LineWidth', 2)
plot([2e3,2e3], ylim, 'r--', 'LineWidth', 2)
xlim([0, 3500])
xlabel('ganch (pN s/mm)')
ylabel('Ring tension (pN)')
xlim([0, 3500])
% legend({'Simulation node velocity', 'Analytical node velocity', ...
%     'Simulation tension', 'Analytical tension'}, 'EdgeColor', 'w', ...
%     'Location', 'NW')