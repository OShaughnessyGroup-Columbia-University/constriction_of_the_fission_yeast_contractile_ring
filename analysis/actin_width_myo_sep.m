% change from original: also do sum intensity of Myo2 and Myp2
% change from v2: use gaussian instead of flat disk
% v3 and v4 skipped
% figure
% change from v5: add calculation of sep_myo_myp and bound_frac
% change from v6: calculate 80% interpercentile of actin 
function [act_spread,sep_myo_myp] = actin_width_myo_sep(rbead,rmyo,rmyp)
%% project components onto z-rho plane
[~,rho_myp] = cart2pol(rmyp(1,:),rmyp(2,:));
z_myp = rmyp(3,:);
[~,rho_myo] = cart2pol(rmyo(1,:),rmyo(2,:));
z_myo = rmyo(3,:);
[~,rho_act] = cart2pol(rbead(1,:),rbead(2,:));
z_act = rbead(3,:);

%% count what fraction of actin is bundled
    rho_myo = median(rho_myo);
    rho_myp = median(rho_myp);
    sep_myo_myp = rho_myo - rho_myp;
    temp = [rbead(3,:);rho_act] - [median(rmyo(3,:));rho_myo];
    temp = sqrt(sum(temp.*temp));
    act_spread = prctile(rho_act,95) - prctile(rho_act,15);