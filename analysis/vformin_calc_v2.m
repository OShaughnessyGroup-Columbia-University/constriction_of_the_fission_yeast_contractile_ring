load gm2_scan_data.mat
load('gm_2/gm_2_1_00sec.mat')
for irun = 1:25
    load(['gm_2/gm_2_' num2str(irun) '_' num2str(600) 'sec.mat'])
    phiform = cart2pol(rbead(1, ifor), rbead(2, ifor));
    phihat = [-sin(phiform) ; cos(phiform); zeros(1,length(phiform))];
    formin_speed = sqrt(sum(vbead(:,ifor).*vbead(:,ifor)));
    for if_idx = 1:length(ifor)
        vphi(irun, if_idx) = phihat(:, if_idx)'*vbead(:, ifor(if_idx));
    end
    gamma_anch(irun) = gm;
end

clear node_speed_cell
for i = 1:5
    node_speed_cell{i} = abs(vphi(i:5:25, :)); 
    vnode_mean(i) = mean(mean(node_speed_cell{i}, 2));
    vnode_std(i) = std(mean(node_speed_cell{i}, 2));
%     vnode_std(i) = std(node_speed_cell{i}(:));

%     tens_ave(i) = mean(tens_store(i:5:end, 21));
%     tens_std(i) = std(tens_store(i:5:end, 21));
end

gamma_anch = gamma_anch(1:5);