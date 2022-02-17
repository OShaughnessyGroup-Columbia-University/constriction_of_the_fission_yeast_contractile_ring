% real_t_min_array = repmat(0,1,21);
real_t_min_array = 0:16;
r_ring_array = 1.85 - 0.07*real_t_min_array;

n_for_exp = [];
n_myo_exp = [];
n_myp_exp = [];
n_for = [];
n_myo = [];
n_myp = [];

load_prefix = 'sizep_17_';
filename = [load_prefix,'00min.mat'];
load(filename);

for icomp = 1:numel(real_t_min_array)
    real_t_min = real_t_min_array(icomp);
    r_ring = r_ring_array(icomp);
    [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev] = modify_binding_rates_v2(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
    n_for_exp = [n_for_exp, 2*pi*r_ring * binding_rate_for / koffmyo];
    n_myo_exp = [n_myo_exp, 2*pi*r_ring * binding_rate_myo2 / koffmyo];
    n_myp_exp = [n_myp_exp, 2*pi*r_ring * binding_rate_myp2 / koffmyp];
    
    filename = [load_prefix,num2str(icomp-1),'min.mat'];
    
    if exist(filename,'file')
        load(filename,'ifor','rmyo','rmyp','r')
    end
    n_for = [n_for, numel(ifor)];
    n_myo = [n_myo, size(rmyo,2)];
    n_myp = [n_myp, size(rmyp,2)];
end

figure
plot(real_t_min_array,n_for_exp,'b-')
hold on
plot(real_t_min_array,n_myo_exp,'r-')
plot(real_t_min_array,n_myp_exp,'g-')
plot(real_t_min_array,n_for,'bo')
plot(real_t_min_array,n_myo,'ro')
plot(real_t_min_array,n_myp,'go')
hold off