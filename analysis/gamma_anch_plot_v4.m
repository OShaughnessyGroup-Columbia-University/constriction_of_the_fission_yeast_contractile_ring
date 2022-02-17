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

gm_theory = 0:50:3500;
% for i = 1:71
%     vnode_theory(i)= mean(lfil.*(rho_myo_list*fhead + rho_myp_list*fheadmyp).*...
%         (1 - (vpol_list - vcont)/vmyo0)./...
%         ((gt(i) + gmsol)*lfil.*rho_myo_list/16 +...
%         lfil.*(2*rho_myo_list*gamma_myo + rho_myp_list*gamma_myp)));
%     t_theory(i) = mean(16*lfil./2.*(rho_myo_list*fhead + rho_myp_list*fheadmyp).*...
%         (1 - (vpol_list-vcont+vnode_theory(i))/vmyo0));
% end

vg_fit_type= fittype('f/(x+g0)', 'independent', 'x');
vg_fit = fit(gamma_anch.', vnode_mean.', vg_fit_type, 'StartPoint', [80, 2000])
vnode_theory = vg_fit(gm_theory);
tv_fit = fit(vnode_mean', tens_ave_lap', 'poly1')
t_theory = tv_fit(vnode_theory);

figure
yyaxis left
set(gca, 'Position',[0.2 0.2 0.6 0.7]);
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
box on
hold on
ylim([0, 60])
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
errorbar(gamma_anch/2e3, tens_ave_lap, tens_std_lap, 'o', 'LineWidth', 2, ...
   'MarkerSize', 12, 'CapSize', 15)
plot(gm_theory/2e3, t_theory, 'LineWidth', 3)
plot(xlim, [640 640], 'LineWidth', 2)
% plot([2e3,2e3], ylim, 'k--', 'LineWidth', 2)
xlim([0, 3500])
xlabel({'Relative node anchor' 'drag coefficient'})
ylabel('Ring tension (pN)')
xlim([0, 1.6])
ylim([0, 800])
% legend({'Simulation node velocity', 'Analytical node velocity', ...
%     'Simulation tension', 'Analytical tension'}, 'EdgeColor', 'w', ...
%     'Location', 'NW')