load('gm_2/gm_2_1_00sec.mat')
% tlist = 0:30:1140;
tlist = 600;
vlist = [];
vphi = nan(25, numel(tlist));
for irun = 4:5:25
% for irun = 1:25
    for istep = 1:length(tlist)
        load(['gm_2/gm_2_' num2str(irun) '_' num2str(tlist(istep)) 'sec.mat'])
        formyo = sum(fmmat);
        phimyo = cart2pol(rmyo(1, :), rmyo(2, :));
        phihat = [-sin(phimyo) ; cos(phimyo); zeros(1,length(phimyo))];
        myo_speed = sqrt( sum( vmyo.*vmyo ) );
        vphi(irun, istep) = mean( diag( phihat'*vmyo ) );
        
        vmlist = diag( phihat'*vmyo );
        histogram(vmlist(formyo == 1))
        
        vlist = [vlist; diag( phihat'*vmyo )];
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