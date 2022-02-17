load('gm_2/gm_2_1_00sec.mat')
for irun = 1:25
    load(['gm_2/gm_2_' num2str(irun) '_600sec.mat'])
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
    vnode_mean(i) = mean(node_speed_cell{i}(:));
    vnode_std(i) = std(node_speed_cell{i}(:));
end

gamma_anch = gamma_anch(1:5);