load('gm_2/gm_2_1_00sec.mat')
% tlist = 0:30:1140;
tlist = 600;
vlist = [];
vphi = nan(25, numel(tlist));
for irun = 1
% for irun = 1:25
    for istep = 1:length(tlist)
        load(['gm_2/gm_2_' num2str(irun) '_' num2str(tlist(istep)) 'sec.mat'])
        phiform = cart2pol(rbead(1, ifor), rbead(2, ifor));
        phihat = [-sin(phiform) ; cos(phiform); zeros(1,length(phiform))];
        formin_speed = sqrt(sum(vbead(:,ifor).*vbead(:,ifor)));
        vphi(irun, istep) = mean(diag(abs(phihat'*vbead(:, ifor))));
        vpol_list(istep) = vpol;
        rho_myo_list(istep) = size(rmyo,2)/(2*pi*r_ring);
        rho_myp_list(istep) = size(rmyp,2)/(2*pi*r_ring);
        lfil(istep) = mean(diff(ifor)-1)*0.1;
        
        numel(diag(phihat'*vbead(:, ifor)))
        vlist = [vlist; diag(phihat'*vbead(:, ifor))];
%         vlist(:, irun) = diag(phihat'*vbead(:, ifor));
    end
    gamma_anch(irun) = gm;
end
node_speed = mean(vphi, 2);


%% load means
for i = 1:5
    vnode_mean(i) = mean(node_speed(i:5:end));
    vnode_std(i) = std(node_speed(i:5:end));
%     vnode_std(i) = std(node_speed_cell{i}(:));

    tens_ave(i) = nanmean(tens{i}(:));
    tens_std(i) = std(nanmean(tens{i}, 2));
    tens_ave_lap(i) = nanmean(tens_lap{i}(:));
    tens_std_lap(i) = std(nanmean(tens_lap{i}, 2));
end

gamma_anch = gamma_anch(1:5);