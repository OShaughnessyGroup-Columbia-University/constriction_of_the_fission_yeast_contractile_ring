%% change from v5: use rk1_v5
% chagne from v6: use rk1_v6
% change from v7: use rk1_v7
% change from v8: use rk1_v8 and parameters_v2
% change from v9: use parameters_v3 and rk1_v9
% change from v10: use rk1_v10 and parameters_v4
% change from v11: use rk1_v11 and parameters_v5
% change from v12: use rk1_v12 and parameters_v6
% v14-17 skipped
% change from v13: use parameters_v11
% change from v18: use parameters_v12
% change from v19: use parameters_v13 and name final3d
% change from v20: use parameters_v14 and name fb2d
% change from v21: use parameters_v15 and name fbp2d
% v23 skipped
% change from v22: use parameters_v16 (no scan), run_pos_z_cyl_v1, constriction, name fcold
% v25 skipped
% change from v24: use rk1_v13 and parameters_v18, name fcnew
% change from v26: use parameters_v29 and name anc2d, output velocities
% change from v41: use parameters_v31 and name oanc1d
% change from v44: use parameters_v35, run_pos_z_cyl, at 12 min, name lnode 
% v45 to v48 skipped
% change from v49: use run_pos_z_cyl_v5, parameters_v36 and name bind2d
% change from v50: parameters_v39, name bind2d4
% v51 to v55 skipped
% change from v56: parameters_v40, name drag2d
% change from v57: parameters_v42, name drag2d2, nstep and dt depend on
% itask
% v58 skipped
% change from v59: use parameters_v43 (same as the one task_yeti_v60 is
% using) and run_pos_z_cyl_v4. full constriction. name drag1d
% v60-62 skipped
% change from v63: use parameters_v47. name drag1d2
% v64 skipped

function [] = task_yeti_v19(itask)
% % make sure that the random function generator can give a different seed
% rng(itask)
% rng(randi(30000))
for i = 1:itask
    a = rand;
end
geometry = 'cyl'
if geometry == 'cyl'
    run_pos_z_cyl_v4
elseif geometry == 'sph'
    run_pos_z_sph
elseif geometry == 'tes'
    test
    geometry = 'cyl';
end
parameters_v47
%% decide which file to load
a = dir;
% first part of the name of the file to load
prefix = 'drag1d2_';
suffix = 'min.mat';

for j = 0:1000
    for i = 1:size(a,1)
        if strcmp(a(i).name, strcat(prefix,num2str(itask),'_',num2str(j),suffix))
            found = true;
            break
        else
            found = false;
        end 
    end
    
    if ~found
        if j == 0
            save(strcat(prefix,num2str(itask),'_0',suffix))
        else
            loadj = j-1;
            break
        end
    end
end

clear j

%% run simulation
totaltime = 40; % ring radius divided by septation rate is 26 minutes  

for t = loadj:loadj + totaltime-1
    load(strcat(prefix,num2str(itask),'_',num2str(t),suffix))
    temp = 100:-1:1;
    nstep = temp(itask)*.5*60; % number of steps per minute
    dt = 60/nstep;
    equil_time_min = 10;
    real_t_min = t * dt * nstep / 60 - equil_time_min; % real time, in minutes

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev] = modify_binding_rates(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
        vs_abs = .07/60;
    end
    
    rk1_v13
    tsave = t+1;
    clear t
    save(strcat(prefix,num2str(itask),'_',num2str(tsave),suffix), 'real_t_min','rbead','rmyo','rmyp','ifor','ipt','bancf','bancm','xmat','fmmat','dbead_first','r','r_ring','fc','fanc','binding_rate_myo2','binding_rate_myp2','binding_rate_for','vpol','rsev','vbead','vmyo')
    if r_ring < 0
        break
    end
end

% load('test_3_1min.mat')
% ringplot_pick
end