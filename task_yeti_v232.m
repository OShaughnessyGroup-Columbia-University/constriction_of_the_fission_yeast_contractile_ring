function [] = task_yeti_v19(itask)
% % make sure that the random function generator can give a different seed
rng(itask)
rng(randi(21000))
for i = 1:itask
    a = rand;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial time
t=0;
debug_cnt = -10000;
texp = 10;

real_t_min = texp;
r = 1.85 - 0.07*real_t_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize simulation
a = dir;
% first part of the name of the file to load
prefix = 'kc2d_';
prefix_late = 'kcrec_'; % for after unbinding threshold is restored
suffix = 'sec.mat';
if(exist(strcat(prefix, num2str(itask), '_00', suffix), 'file'))
    load(strcat(prefix, num2str(itask), '_00', suffix))
else
    geometry = 'cyl';
    run_pos_z_cyl_v8
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
parameters_v95

x1 = [40,  5, 40];
x2 = [ 5, 30, 30];
[X1, ~] = meshgrid(x1, 1:10);
[X2, ~] = meshgrid(x2, 1:10);

myo2p_unbinding_force = X1(itask);
myp2p_unbinding_force = X2(itask);
kcap = myo2p_unbinding_force / (rcapmyo_short-r_flat_myo);
kcap_p = myp2p_unbinding_force / (rcapmyp-r_flat_myp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look for files to load
tprint = 20;
totaltime = 10*60; % ring radius divided by septation rate is 26 minutes  
for j = 0:tprint:totaltime
    disp(strcat(prefix, num2str(itask), '_', num2str(j), suffix))
    if exist(strcat(prefix, num2str(itask), '_', num2str(j), suffix), 'file')
        loadj = j
	found = true
    elseif exist(strcat(prefix_late, num2str(itask), '_', num2str(j), suffix),'file')
        prefix = prefix_late
        loadj = j
	found = true
    else
        found = false
    end
    
    if ~found
        if j == 0
            save(strcat(prefix, num2str(itask), '_0', suffix))
        else
            loadj = j-tprint;
            break
        end
    end
end

clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run simulation
totaltime = 10*60; % ring radius divided by septation rate is 26 minutes  
warmuptime = .5;   % in seconds
equiltime = -10;   % in seconds
t=-equiltime;

debug_cnt = -10000;
if loadj==0
    disp('Warming up')
    load(strcat(prefix, num2str(itask), '_0', suffix))
    save(strcat(prefix, num2str(itask), '_00', suffix))
    t_rbead = rand(1, length(rbead))*1e-6;
    t_rmyo = rand(1, length(rmyo))*1e-6;
    t_rmyp = rand(1, length(rmyp))*1e-6;
    dt = .005; 
    nstep = warmuptime/dt;
    real_t_min = texp;

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev]=modify_binding_rates_v2(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
        vs_abs = 0;
    end
    rk1_v48gauss

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Equilibrating')
    for et = 0:tprint:equiltime
        disp(['Equilibration: t=' num2str(et) ' seconds'])
        t = et - equiltime;
        dt = 1./30; 
        nstep = tprint/dt;
        rk1_v48gauss
	save(strcat(prefix, num2str(itask), '_', num2str(et), 'e', suffix), 'real_t_min', 'rbead', 'rmyo', ...
			'rmyp', 'bmmat', 'bpmat', 'ifor', 'ipt', 'bancf', 'bancm', 'xmat', 'fmmat', 'dbead_first', ...
			'r', 'r_ring', 'fc', 'fanc', 'binding_rate_myo2', 'binding_rate_myp2', 'binding_rate_for', ...
			'vpol', 'rsev', 'vbead', 'vmyo', 'vmyp', 't_rbead', 't_rmyo', 't_rmyp');
    end
    save(strcat(prefix, num2str(itask), '_0', suffix), 'real_t_min', 'rbead', 'rmyo', ...
        'rmyp', 'bmmat', 'bpmat', 'ifor', 'ipt', 'bancf', 'bancm', 'xmat', 'fmmat', 'dbead_first', ...
        'r', 'r_ring', 'fc', 'fanc', 'binding_rate_myo2', 'binding_rate_myp2', 'binding_rate_for', ...
        'vpol', 'rsev', 'vbead', 'vmyo', 'vmyp', 't_rbead', 't_rmyo', 't_rmyp');
end

debug_cnt = -1;
if debug_cnt == 0
	status = mkdir('dump');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Commencing run')
for t = loadj:tprint:loadj + totaltime-1
    disp(['t = ' num2str(t)])
    load(strcat(prefix, num2str(itask), '_', num2str(t), suffix))
    if(t >= 180)
        prefix = prefix_late;
        myo2p_unbinding_force = 40;
        kcap = myo2p_unbinding_force / (rcapmyo_short-r_flat_myo);
        myp2p_unbinding_force = 30;
        kcap_p = myp2p_unbinding_force / (rcapmyp-r_flat_myp);
        vs_abs = 0.07/60;
    else
        vs_abs = 0;
    end 
    dt = 1/30.;
    nstep = tprint/dt;
    real_t_min = texp;

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev]=modify_binding_rates_v2(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
    end
    rk1_v48gauss
    tsave = t+tprint;
    clear t
    save(strcat(prefix, num2str(itask), '_', num2str(tsave), suffix), 'real_t_min', 'rbead', ...
        'rmyo', 'rmyp', 'bmmat', 'bpmat', 'ifor', 'ipt', 'bancf', 'bancm', 'xmat', 'fmmat', ...
        'dbead_first', 'r', 'r_ring', 'fc', 'fanc', 'binding_rate_myo2', 'binding_rate_myp2', ...
        'binding_rate_for', 'vpol', 'rsev', 'vbead', 'vmyo', 'vmyp', 't_rbead', 't_rmyo', 't_rmyp');
    if r_ring < 0
        break
    end
end
end
