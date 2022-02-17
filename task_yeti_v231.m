function [] = task_yeti_v19(itask)
% % make sure that the random function generator can give a different seed
rng(itask)
rng(randi(21000))
for i = 1:itask
    a = rand;
end
geometry = 'cyl'
tfix = 0;
real_t_min = tfix;
r = 1.85 - real_t_min*0.07;
if geometry == 'cyl'
    run_pos_z_cyl_v8
elseif geometry == 'sph'
    run_pos_z_sph
elseif geometry == 'tes'
    test
    geometry = 'cyl';
end
parameters_v95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 1./[600, 100, 30, 10];
x2 = [0, 10, 20, 40, 100, 200];
x3 = 1:10;
[X1, X2, ~] = meshgrid(x1, x2, x3);
dtrun = X1(itask);
ga = X2(itask);
ga_myo2 = ga;
ga_myp2 = ga;

%% decide which file to load
a = dir;
% first part of the name of the file to load
prefix = 'gaxdt_'
suffix = 'sec.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tprint = 10;
totaltime = 3*60; % ring radius divided by septation rate is 26 minutes  
for j = 0:tprint:totaltime
    for i = 1:size(a, 1)
        if strcmp(a(i).name, strcat(prefix, num2str(itask), '_', num2str(j), suffix))
            loadj = j;
            found = true;
            break
        else
            found = false;
        end 
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

%% run simulation
totaltime = 3*60; % ring radius divided by septation rate is 26 minutes  
warmuptime = .5;   % in seconds
equiltime = -1;   % in seconds
t=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug_cnt = -10000;
if loadj==0
    disp('Warming up')
    load(strcat(prefix, num2str(itask), '_0', suffix))
    save(strcat(prefix, num2str(itask), '_00', suffix))
    t_rbead = rand(1, length(rbead))*1e-6;
    t_rmyo = rand(1, length(rmyo))*1e-6;
    t_rmyp = rand(1, length(rmyp))*1e-6;
    dt = dtrun;
    nstep = warmuptime/dt;
    real_t_min = tfix;

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
        dt = dtrun;
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
    dt = dtrun
    nstep = tprint/dt;
    real_t_min = tfix;
    %real_t_min = t/60;

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev]=modify_binding_rates_v2(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
        vs_abs = 0;%0.07/60;
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
