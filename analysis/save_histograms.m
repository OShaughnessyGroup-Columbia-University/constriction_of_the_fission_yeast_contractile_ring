%% prepare data structures
% clear vm_med vp_med vm_sd vp_sd mp_frac pm_frac mm_frac pp_frac vmod
pfx = 'gm5';

load([pfx '_1_00sec.mat']);
nval2 = 1 %size(X1, 2); % values of x1
nval1 = size(X1, 1); % values of x2
nval = nval1*nval2

ilist = 5;
% ilist = 1:nval:11*nval;
irmax = 20*nval;
tlist = 40:10:240;
v_store = nan(irmax, numel(tlist));
tf_store = nan(irmax, numel(tlist));
ts_store = nan(irmax, numel(tlist));
tlap_store = nan(irmax, numel(tlist));
tcirc4_store = nan(irmax, numel(tlist));

%% calculate and save histogramsw
for i0 = ilist
    run_list = i0:nval:irmax;
%     run_list = i0:i0+nval-1;
    vagg = [];
    pagg = [];
    disp(['* * * * * * * * * * * * i0 = ' num2str(i0) ' * * * * * * * * * * * *'])
    for it = 1:length(tlist)
        trun = tlist(it);
        disp(['t = ' num2str(trun) ' sec'])
        %% aggregate data and calculate tension for this drag value
        for idx = 1:length(run_list)
            irun = run_list(idx);
            flag = 0;
%             disp(['irun = ' num2str(irun)])
            try
                load([pfx '_' num2str(irun) '_' num2str(trun) 'sec.mat']);
                flag = 1;
%                 vplist(i0) = vpol;
%                 gmlist(i0) = gm;
            catch
%                 disp(['No file waxv4_' num2str(irun) '_' num2str(trun) 'sec.mat'])
                continue
            end
            
            if(flag == 1)
                %% calculating tension contributions and velocity
                tension_laplace_v1
                tlap_store(irun, it) = mean(tens);
                tension_circ_v4
                tcirc4_store(irun, it) = mean(tens);
                calc_sliding_fixed

                tf_store(irun, it) = mean(tfixed);
                ts_store(irun, it) = mean(tslide);
                v_store(irun, it) = vave;

                % saving polarity and velocity data
                pagg = [pagg p];
                vagg = [vagg v];
            end
        end
    end
    
    if(isempty(pagg))
        disp('No runs found with these parameters!')
        continue
    end
%     save(['./histograms/gm' num2str(i0) '_waxv4e.mat'])
    %% save histograms
    bins = -100:4:100;
    figure
    hold on
    xlabel('Node velocity (nm/s)')
    ylabel('Probability')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
%     title(['v_{pol} = ' num2str(1e3*vplist(i0)) ' nm/s'])
    title([pfx ', i0=' num2str(i0)])
%     mnodes = pagg(1, :) < pagg(2,:);
%     pnodes = pagg(1, :) > pagg(2,:);
%     fnodes = pagg(1, :) ~= 0 | pagg(2,:) ~= 0;
    mnodes = pagg(1, :) == 0 & pagg(2,:) >= 1;
    pnodes = pagg(1, :) >= 1 & pagg(2,:) == 0;
    fnodes = pagg(1, :) ~= 0 | pagg(2,:) ~= 0;
%     histogram(1e3*vagg(pnodes), bins, 'Displayname', 'p>0', 'normalization', 'probability')
%     histogram(1e3*vagg(mnodes), bins, 'Displayname', 'p<0', 'normalization', 'probability')
%     for p1 = 0:4
%         for p2 = 0:4
%             histogram(1e3*vagg(pagg(1, :) == p1 && pagg(2, :) == p2), bins, ...
%                 'Displayname', ['p=' num2str(p1) ',' num2str(p2)], 'normalization', 'probability')
%         end
%     end
    
    for nf = 0:4
        histogram(1e3*vagg(sum(pagg) == nf), bins, ...
            'Displayname', ['nf=' num2str(nf)], 'normalization', 'probability')  
    end
    
    histogram(1e3*vagg(sum(pagg) ~= 0), bins, 'Displayname', 'All nodes', 'normalization', 'probability')
    cts = zeros*bins(2:end);
    for nfor = 1:4
        cts = cts + nfor*histcounts(1e3*vagg(sum(pagg)==nfor), bins);
    end
%     bar((bins(2:end)+bins(1:end-1))/2, cts/sum(cts), 'Displayname', 'All forsec')
%     histogram(1e3*vagg(pagg(1, :) == pagg(2,:)), bins, 'Displayname', 'p=0', 'normalization', 'probability')
%     disp('average mod velocity of forsecs')
    vmod(i0) = mean(1e3*abs(vagg(mnodes | pnodes)));
    vmod_sd(i0) = std(1e3*abs(vagg(mnodes | pnodes)));
    disp(['mean velocity of (-) nodes ' num2str(mean(1e3*vagg(mnodes)))])
    vm_med(i0) = mean(1e3*vagg(mnodes));
    vm_sd(i0) = std(1e3*vagg(mnodes));
    disp(['mean velocity of (+) nodes ' num2str(mean(1e3*vagg(pnodes)))])
    vp_med(i0) = mean(1e3*vagg(pnodes));
    vp_sd(i0) = std(1e3*vagg(pnodes));
    sarle(i0) = (vp_med(i0) - vm_med(i0))/vp_sd(i0);
%     sarle(i0) = (1+skewness(pnodes | mnodes)^2)/kurtosis(pnodes| mnodes);
    disp(['bimodality of vfor distribution:' num2str(sarle(i0))])
%     disp('median velocity of null nodes')
%     v0_med(i0) = mean(pagg(1, :) == 1 & pagg(2,:) == 0))
%     disp('fraction of (-) nodes going in (-) dirn:')
    mm_frac(i0) = sum(vagg(mnodes) < 0)/length(vagg(mnodes));
%     disp('fraction of (-) nodes going in (+) dirn:')
    mp_frac(i0) = sum(vagg(mnodes) > 0)/length(vagg(mnodes));
%     disp('fraction of (+) nodes going in (-) dirn:')
    pm_frac(i0) = sum(vagg(pnodes) < 0)/length(vagg(pnodes));
%     disp('fraction of (+) nodes going in (+) dirn:')
    pp_frac(i0) = sum(vagg(pnodes) > 0)/length(vagg(pnodes));
    saveas(gcf, ['plots/' pfx ', i0=' num2str(i0) '.png']);
    savefig(gcf, ['plots/' pfx ', i0=' num2str(i0) '.fig']);
%     save([pfx '_i0' num2str(i0) '.mat'])
%     close
%     save(['waxv4_vp' num2str(vplist(i0)) '.mat'])
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

% bins = -50:2:50;
% figure
% hold on
% set(gca, 'Linewidth', 3)
% set(gca, 'FontSize', 24)
% set(gca, 'Position', [0.2 0.2 0.7 0.7]);
% for i = 0:1
%     for j = 0:1
%         if i+j <=2
%             %title(['(' num2str(i) ',' num2str(j) '), \gamma = ' num2str(i0/10) 'x wt'])
%             xlabel('Node velocity (nm/s)')
%             ylabel('Probability')
%             %ylim([0 0.25])
%             histogram(1e3*vagg(sum(pagg == [i; j]) == 2), bins, 'Displayname', ['n(+)=' num2str(i) ', n(-)=' num2str(j)], 'normalization', 'probability')
%         end
%     end
% end