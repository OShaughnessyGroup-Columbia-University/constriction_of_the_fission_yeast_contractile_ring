% change from the original: binding_rate_myo2 and binding_rate_myp2 match binding_rate_myo2_init
% and binding_rate_myp2_init at real_t_min = 0
% change from v1: smooth formin curve
function [binding_rate_myo2, binding_rate_myp2, binding_rate_cdc12, vpol, rsev] = modify_binding_rates(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo)
%     if real_t_min < 5
%         binding_rate_myo2 = binding_rate_myo2_init * (183 + 9/5 * real_t_min)/183;
%         binding_rate_myp2 = binding_rate_myp2_init * (241 + 48/5 * real_t_min)/241;
%     elseif real_t_min < 15
%         binding_rate_myo2 = binding_rate_myo2_init * (192 + 84/10 * (real_t_min-5))/183;
%         binding_rate_myp2 = binding_rate_myp2_init * (289 + 270/10 * (real_t_min-5))/241;
%     else
%         binding_rate_myo2 = binding_rate_myo2_init * (276 + 181/10 * (real_t_min-15))/183;
%         binding_rate_myp2 = binding_rate_myp2_init * (559 + 41/10 * (real_t_min-15))/241;
%     end
    
    time_wu = [-5, 5, 15, 25];
    rho_myo2_exp = [18, 22, 32, 52];
    rho_myp2_exp = [9, 17, 32, 34];
    binding_rate_myo2 = spline(time_wu, rho_myo2_exp, real_t_min)/spline(time_wu, rho_myo2_exp, 0) * binding_rate_myo2_init;
    binding_rate_myp2 = spline(time_wu, rho_myp2_exp, real_t_min)/spline(time_wu, rho_myp2_exp, 0) * binding_rate_myp2_init;
    
    time_chen = 0:2:26;
    n_cdc12_exp = [400, 383, 383, 383, 373, 363, 376, 389, 358, 347, 319, 300, 248, 186];
    n_for3_exp = [300, 311, 323, 322, 322, 286, 305, 285, 277, 212, 185, 173, 150, 127];
    temp = smooth(spline(time_chen, n_cdc12_exp, 0:26));
    n_cdc12 = temp(round(real_t_min)+1);
    binding_rate_cdc12 = n_cdc12/400 * (r/r_ring) * binding_rate_for_init;
    
    rsev_vpol_fit_noplot
end