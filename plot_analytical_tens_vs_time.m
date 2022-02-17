% real_t_min_array = repmat(0,1,21);
real_t_min_array = 0:16;
r_ring_array = 1.85 - 0.07*real_t_min_array;

n_for_exp = [];
n_myo_exp = [];
n_myp_exp = [];
vpol_exp = [];
l_exp = [];

for icomp = 1:numel(real_t_min_array)
    real_t_min = real_t_min_array(icomp);
    r_ring = r_ring_array(icomp);
    [binding_rate_myo2, binding_rate_myp2, binding_rate_cdc12, vpol, rsev, l_mean_goal] = modify_binding_rates_v3(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
    n_for_exp = [n_for_exp, 2*pi*r_ring * binding_rate_cdc12 / koffmyo];
    n_myo_exp = [n_myo_exp, 2*pi*r_ring * binding_rate_myo2 / koffmyo];
    n_myp_exp = [n_myp_exp, 2*pi*r_ring * binding_rate_myp2 / koffmyp];
    vpol_exp = [vpol_exp, vpol];
    l_exp = [l_exp, l_mean_goal];
%     filename = [load_prefix,num2str(icomp-1),'min.mat'];
%     
%     if exist(filename,'file')
%         load(filename,'ifor','rmyo','rmyp','r')
%     end
%     n_for = [n_for, numel(ifor)];
%     n_myo = [n_myo, size(rmyo,2)];
%     n_myp = [n_myp, size(rmyp,2)];
end

% figure
% plot(real_t_min_array,n_for_exp,'b-')
% hold on
% plot(real_t_min_array,n_myo_exp,'r-')
% plot(real_t_min_array,n_myp_exp,'g-')
% % plot(real_t_min_array,n_for,'bo')
% % plot(real_t_min_array,n_myo,'ro')
% % plot(real_t_min_array,n_myp,'go')
% hold off

tens_analytical = 16*1.3 * (1-vpol_exp/.24) .* (n_myo_exp+n_myp_exp)./r_ring_array/2/pi .* l_exp/2;
phi = acosd(r_ring_array / 1.85);
plot(phi,tens_analytical,'r-')
axis([0,80,0,1500])
xlabel('Ring position, \phi (deg)')
ylabel('Ring tension (pN)')