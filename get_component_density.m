% Calculate experimental densities of various components relative to onset time
function [rhom, rhop, rhof, vpol, rsev] = get_component_density(real_t_min, r_ring, koffmyo)
    time_wu = [-5, 5, 15, 25];
    rho_myo2_exp = [18, 22, 32, 52];
    rho_myp2_exp = [9, 17, 32, 34];
    rhom = spline(time_wu, rho_myo2_exp, real_t_min)/spline(time_wu, rho_myo2_exp, 0);
    rhop = spline(time_wu, rho_myp2_exp, real_t_min)/spline(time_wu, rho_myp2_exp, 0);
    
    time_chen = 0:2:26;
    n_cdc12_exp = [400, 383, 383, 383, 373, 363, 376, 389, 358, 347, 319, 300, 248, 186];
    temp = smooth(spline(time_chen, n_cdc12_exp, 0:26));
    n_cdc12 = temp(round(real_t_min)+1);
    rhof = n_cdc12/r_ring;
    
    rsev_vpol_fit_noplot
end
