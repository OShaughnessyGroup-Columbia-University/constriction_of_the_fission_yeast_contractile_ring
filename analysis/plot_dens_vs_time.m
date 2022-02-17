% real_t_min_array = repmat(0,1,21);
real_t_min_array = 0:16;
r_ring_array = 1.85 - 0.07*real_t_min_array;

rho_for_exp = [];
rho_myo_exp = [];
rho_myp_exp = [];
rho_for = [];
rho_myo = [];
rho_myp = [];

load_prefix = 'fh2d5_50_';
filename = [load_prefix,'00min.mat'];
load(filename);

for icomp = 1:numel(real_t_min_array)
    real_t_min = real_t_min_array(icomp);
    r_ring = r_ring_array(icomp);
    [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev] = modify_binding_rates_v2(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
    rho_for_exp = [rho_for_exp, binding_rate_for / koffmyo];
    rho_myo_exp = [rho_myo_exp, binding_rate_myo2 / koffmyo];
    rho_myp_exp = [rho_myp_exp, binding_rate_myp2 / koffmyp];
    
    filename = [load_prefix,num2str(icomp-1),'min.mat'];
    
    if exist(filename,'file')
        load(filename,'ifor','rmyo','rmyp','r')
    end
    rho_for = [rho_for, numel(ifor)/2/pi/r_ring];
    rho_myo = [rho_myo, size(rmyo,2)/2/pi/r_ring];
    rho_myp = [rho_myp, size(rmyp,2)/2/pi/r_ring];
end

figure
plot(real_t_min_array,rho_for_exp,'b-')
hold on
plot(real_t_min_array,rho_myo_exp,'r-')
plot(real_t_min_array,rho_myp_exp,'g-')
plot(real_t_min_array,rho_for,'bo')
plot(real_t_min_array,rho_myo,'ro')
plot(real_t_min_array,rho_myp,'go')
hold off