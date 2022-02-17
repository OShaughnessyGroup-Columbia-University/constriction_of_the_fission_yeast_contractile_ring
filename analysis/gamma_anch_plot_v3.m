% figure
% load('gm_2/gm_2_10_600sec.mat')
tens_sem = tens_std/sqrt(5*41-1);
vpol_mean = mean(vpol_list);
rho_myo = mean(rho_myo_list);
rho_myp = mean(rho_myp_list);
l = mean(lfil);

gm_theory = 0:100:3500;
gamma_myo = fhead/vmyo0;
gamma_myp = fheadmyp/vmyo0;
vcont = 2*pi*0.07/60;
% vnode_theory = l*(rho_myo*fhead + rho_myp*fheadmyp)*(1 - (vpol_mean - vcont)/vmyo0)./...
%     ((gm_theory + gmsol) + l*(2*rho_myo*gamma_myo + rho_myp*gamma_myp));
vnode_theory = l*(rho_myo*fhead + rho_myp*fheadmyp)*(1 - (vpol - vcont)/vmyo0)./...
    ((gm_theory + gmsol)*l*rho_myo/16 + l*(2*rho_myo*gamma_myo + rho_myp*gamma_myp));
% vnode_theory = l*(rho_myo*fhead + rho_myp*fheadmyp)*(1 - (vpol - vcont)/vmyo0)./...
%     ((gm_theory/2.5 + gmsol)*l*rho_myo/16 + l*(2*rho_myo*gamma_myo + rho_myp*gamma_myp));
t_theory = 16*l/2*(rho_myo*fhead + rho_myp*fheadmyp)*...
    (1 - (vpol_mean-vcont+vnode_theory)/vmyo0);
% tens_analytical = 16*1.3 * (1-vpol_exp/.24) .* (n_myo_exp+n_myp_exp)...
%     ./r_ring_array/2/pi .* l_exp/2;

figure
yyaxis left
set(gca, 'Position',[0.2 0.2 0.6 0.7]);
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
box on
hold on
ylim([0, 80])
errorbar(gamma_anch/2e3, 1e3*vnode_mean, 1e3*vnode_std, 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'auto', 'MarkerSize', 12, 'CapSize', 15)
plot(gm_theory/2e3, 1e3*vnode_theory, 'LineWidth', 3)
plot(xlim, [28, 28], 'LineWidth', 2)
plot([1, 1], ylim, 'k--', 'LineWidth', 2)
ylabel('vnode (nm/s)')
xlabel({'Relative node anchor' 'drag coefficient'})
xlim([0, 1.6])

yyaxis right
set(gca, 'Position',[0.2 0.2 0.6 0.7]);
set(gca, 'Linewidth', 3)
set(gca, 'FontSize', 30)
box on
hold on
% errorbar(gamma_anch/2e3, tens_ave, tens_sem, 'o', 'LineWidth', 2, ...
%    'MarkerSize', 12, 'CapSize', 15)
errorbar(gamma_anch/2e3, tens_ave, tens_std, 'o', 'LineWidth', 2, ...
   'MarkerSize', 12, 'CapSize', 15)
plot(gm_theory/2e3, t_theory, 'LineWidth', 3)
plot(xlim, [640 640], 'LineWidth', 2)
% plot([2e3,2e3], ylim, 'k--', 'LineWidth', 2)
xlim([0, 3500])
xlabel({'Relative node anchor' 'drag coefficient'})
ylabel('Ring tension (pN)')
xlim([0, 1.6])
ylim([0, 900])
% legend({'Simulation node velocity', 'Analytical node velocity', ...
%     'Simulation tension', 'Analytical tension'}, 'EdgeColor', 'w', ...
%     'Location', 'NW')