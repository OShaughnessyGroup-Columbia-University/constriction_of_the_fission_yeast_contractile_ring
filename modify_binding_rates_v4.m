% change from the original: binding_rate_myo2 and binding_rate_myp2 match binding_rate_myo2_init
% and binding_rate_myp2_init at real_t_min = 0
% change from v1: smooth formin curve
function [binding_rate_myo2, binding_rate_myp2, binding_rate_cdc12, vpol, rsev] = ...
    modify_binding_rates_v4(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init,...
    binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo)
    time_wu = [-20, -17, -12, -5, 5, 15, 25];
    rho_myo2_exp = [18, 18, 18, 18, 22, 32, 52];
    %rho_myo2_exp = [9, 12, 16, 18, 22, 32, 52];
    rho_myp2_exp = [0, 0, 4, 9, 17, 32, 34];
    binding_rate_myo2 = spline(time_wu, rho_myo2_exp, real_t_min)/spline(time_wu, rho_myo2_exp, -20) ...
        * binding_rate_myo2_init;
    binding_rate_myp2 = 0*(real_t_min < -17) + spline(time_wu, rho_myp2_exp, real_t_min)/...
        spline(time_wu, rho_myp2_exp, 0) * binding_rate_myp2_init * spline(time_wu, rho_myo2_exp, 0) / ...
        spline(time_wu, rho_myo2_exp, -20) * (real_t_min >= -17);
    disp(['real t min: ' num2str(real_t_min)])
    disp(['binding_rate_myp2: ' num2str(binding_rate_myp2)])
    
    time_chen = [-20, -10:2:26];
    n_cdc12_exp = [360, 360, 369, 387, 375, 396, 400, 383, 383, 383, 373, 363, 376, 389, 358, 347, 319, 300, 248, 186];
    n_for3_exp = [295, 295, 310, 300, 292, 290, 300, 311, 323, 322, 322, 286, 305, 285, 277, 212, 185, 173, 150, 127];
    temp = smooth(spline(time_chen, n_cdc12_exp, -20:26));
    n_cdc12 = temp(round(real_t_min)+21);
    binding_rate_cdc12 = n_cdc12/360 * (r/r_ring) * binding_rate_for_init;
    
    rsev_vpol_fit_noplot_v2
end
