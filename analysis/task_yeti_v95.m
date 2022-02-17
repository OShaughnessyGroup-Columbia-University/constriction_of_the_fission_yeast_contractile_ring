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
% change from v59: use dt = (.2/3) s, name drag2d3
% v60-68 skipped
% change from v69: use rk1_v18, name drag2d4
% change from v70: use dt = .2 s, name drag2d5
% change from v71: use rk1_v19
% v72 skipped
% change from v73: use warm up, name wmup2
% change from v74: use rk1_v20, name wmup3
% change from v75: use parameters_v51, name wmup4
% change from v76: use rk1_v21, name wmup5
% change from v77: use rk1_v22
% change from v78: use rk1_v23
% change from v79: use rk1_v24
% change from v80: use parameters_v52
% change from v81: use rk1_v25
% change from v82: use parameters_v53
% change from v83: use rk1_v26
% change from v84: use rk1_v27
% change from v85: use rk1_v28, name wmup6
% change from v86: use rk1_v29, name myok
% change from v87: use rk1_v30, name mypuni
% change from v88: use parameters_v54, name nosep
% change from v89: use parameters_v55, name nosep2
% change from v90: use parameters_v57, run_pos_z_cyl_v7, name nosep3
% v91 skipped
% change from v92: use run_pos_z_cyl_v8, full constriction, no initial 10 min, name nsfc
% change from v93: use modify_binding_rates_v1, name nsfc2
% v94 skipped
function [] = task_yeti_v19(itask)
% % make sure that the random function generator can give a different seed
rng(itask)
rng(randi(30000))
for i = 1:itask
    a = rand;
end
geometry = 'cyl'
if geometry == 'cyl'
    run_pos_z_cyl_v8
elseif geometry == 'sph'
    run_pos_z_sph
elseif geometry == 'tes'
    test
    geometry = 'cyl';
end
parameters_v57
%% decide which file to load
a = dir;
% first part of the name of the file to load
prefix = 'nsfc2_';
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
totaltime = 30; % ring radius divided by septation rate is 26 minutes  
warmuptime = .5;

if loadj==0
    load(strcat(prefix,num2str(itask),'_0',suffix))
    dt = .005; 
    nstep = warmuptime/dt;
    real_t_min = 0;

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev] = modify_binding_rates(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
        vs_abs = 0;
    end
    rk1_v30
%     save(strcat(prefix,num2str(itask),'_0',suffix), 'real_t_min','rbead','rmyo','rmyp','ifor','ipt','bancf','bancm','xmat','fmmat','dbead_first','r','r_ring','fc','fanc','binding_rate_myo2','binding_rate_myp2','binding_rate_for','vpol','rsev','vbead','vmyo')
end
for t = loadj:loadj + totaltime-1
    load(strcat(prefix,num2str(itask),'_',num2str(t),suffix))
    dt = .2;
    nstep = 5*60;
    equil_time_min = 0;
    real_t_min = t * dt * nstep / 60 - equil_time_min; % real time, in minutes

    if real_t_min >=0
        [binding_rate_myo2, binding_rate_myp2, binding_rate_for, vpol, rsev] = modify_binding_rates_v1(real_t_min, binding_rate_myo2_init, binding_rate_myp2_init, binding_rate_for_init, kofffor, r, r_ring, rhof, koffmyo);
        % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
        vs_abs = .07/60;
    end
    rk1_v30
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