%% prepare data structures
nval = 20;
ilist = 5:5:nval;
tlist = 60:10:600;
v_store = nan(numel(ilist), numel(tlist));
tf_store = nan(numel(ilist), numel(tlist));
ts_store = nan(numel(ilist), numel(tlist));
tlap_store = nan(numel(ilist), numel(tlist));
tcirc4_store = nan(numel(ilist), numel(tlist));

delt = 10.;
load(['noagm5_21_00sec.mat']);
clear vsmooth
for i0 = ilist    
    vagg = [];
    pagg = [];
    disp(['* * * * * * * * * * * * i0 = ' num2str(i0) ' * * * * * * * * * * * *'])
    for it = 1:length(tlist)
        trun = tlist(it);
%         disp(['t = ' num2str(trun) ' sec'])
        %% aggregate data and calculate tension for this drag value
        for irun = i0:nval:200
%             disp(['irun = ' num2str(irun)])
            %% load snapshot delt seconds into the future
            try
                load(['noagm5_' num2str(irun) '_' num2str(trun+delt) 'sec.mat']);
%                 rmyo2 = rmyo;
                phi2myo = cart2pol(rmyo(1, :), rmyo(2, :));
                t2 = t_rmyo;
                load(['noagm5_' num2str(irun) '_' num2str(trun) 'sec.mat']);
                phimyo = cart2pol(rmyo(1, :), rmyo(2, :));
                
                for inode = 1:length(t_rmyo)
                    idx_new = find(t_rmyo(inode)==t2);
                    if(isempty(idx_new))
%                         vsmooth(:, inode) = nan(3, 1);
                        vsmooth(inode) = nan;
                    else
%                         vsmooth(:, inode) = (rmyo2(:, idx_new) - rmyo(:, inode))/delt;
                        dphi = phi2myo(idx_new) - phimyo(inode);
                        dphi = dphi - 2*pi*(dphi > pi) + 2*pi*(dphi < -pi);
                        vsmooth(inode) = r_ring*dphi/delt;
                    end
                end
                gmlist(i0) = gm;
%                 vmyo = (rmyo2 - rmyo1)/(nstep*dt);
%                 vmyp = (rmyp2 - rmyp1)/(nstep*dt);
%                 vbead = (rbead2 - rbead1)/(nstep*dt);
            catch
%                 disp('No file')
                continue
            end
            
%             vmyo = vsmooth;
            %% calculating tension contributions and velocity
%             tension_laplace_v1
%             tlap_store(irun, it) = mean(tens);
%             tension_circ_v4
%             tcirc4_store(irun, it) = mean(tens);
            calc_sliding_fixed
            
            tf_store(irun, it) = mean(tfixed);
            ts_store(irun, it) = mean(tslide);
            v_store(irun, it) = vave;
            
            % saving polarity and velocity data
            pagg = [pagg p];
            vagg = [vagg vsmooth];
        end
    end
    
%     save(['./histograms/gm' num2str(i0) '_gm_2e.mat'])
    %% save histograms
    bins = -100:2:100;
    figure
    hold on
    xlabel('Node velocity (nm/s)')
    ylabel('Probability')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    title(['\gamma = ' num2str(gmlist(i0)) ' pN s/\mum'])
    mnodes = pagg(1, :) == 0 & pagg(2,:) == 1;
    pnodes = pagg(1, :) == 1 & pagg(2,:) == 0;
    histogram(1e3*vagg(pnodes), bins, 'Displayname', 'p>0', 'normalization', 'probability')
    histogram(1e3*vagg(mnodes), bins, 'Displayname', 'p<0', 'normalization', 'probability')
%     histogram(1e3*vagg(pagg(1, :) == pagg(2,:)), bins, 'Displayname', 'p=0', 'normalization', 'probability')
    disp('mean modulus of velocity')
    vmod(i0) = nanmean(abs(vagg)*1e3)
    disp('median velocity of (-) nodes')
    vm_med(i0) = nanmean(1e3*vagg(mnodes))
    disp('median velocity of (+) nodes')
    vp_med(i0) = nanmean(1e3*vagg(pnodes))
%     disp('median velocity of null nodes')
%     v0_med(i0) = mean(pagg(1, :) == 1 & pagg(2,:) == 0))
    disp('fraction of (-) nodes going in (-) dirn:')
    mm_frac(i0) = sum(vagg(mnodes) < 0)/length(vagg(mnodes))
    disp('fraction of (-) nodes going in (+) dirn:')
    mp_frac(i0) = sum(vagg(mnodes) > 0)/length(vagg(mnodes))
    disp('fraction of (+) nodes going in (-) dirn:')
    pm_frac(i0) = sum(vagg(pnodes) < 0)/length(vagg(pnodes))
    disp('fraction of (+) nodes going in (+) dirn:')
    pm_frac(i0) = sum(vagg(pnodes) > 0)/length(vagg(pnodes))
    save(['noagm5_gm' num2str(gmlist(i0)) '.mat'])
    
%     for i = 0:4
%         for j = 0:4
%             if i+j <=4
%                 figure
%                 hold on
%                 title(['(' num2str(i) ',' num2str(j) '), \gamma = ' num2str(gmlist(i0)) 'x wt'])
%                 xlabel('Node velocity (nm/s)')
%                 ylabel('Probability')
%                 ylim([0 0.25])
%                 histogram(1e3*vagg(sum(pagg == [i; j]) == 2), bins, 'Displayname', ['n(+)=' num2str(i) ', n(-)=' num2str(j)], 'normalization', 'probability')
%                 saveas(gcf, ['./histograms/gm' num2str(i0/10) '_(' num2str(i) ',' num2str(j) ').png'])
%                 close
%             end
%         end
%     end
end