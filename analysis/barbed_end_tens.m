pfx = 'pio8kx';
% load([pfx '_1_00min.mat'])
load([pfx '_1_00sec.mat'])

bdmax = 60-1;
nval = 25;
imax = 10*nval;
it = 90;
pback = zeros(1, nval);
pcomp = zeros(1, nval);
x = (0:bdmax-1)*0.1;
% figure(100)
% hold on
for i0 = 1:nval
% for i0 = good
%     figure(100)
    ilist = i0:nval:imax;
    tb_store = [];
    l_store = [];
    vf_store = [];
    isx_store = [];
    tseg_store = nan(length(ilist), bdmax-1);
    
    for idx = 1:length(ilist)
        irun = ilist(idx);
        try
            load([pfx '_' num2str(irun) '_' num2str(it) 'sec.mat'])
%             load([pfx '_' num2str(irun) '_' num2str(it) 'min.mat'])
        catch
            disp([pfx '_' num2str(irun) '_' num2str(it) 'min.mat does not exist'])
            continue
        end
        xlist = find(sum(xmat));
        isx = zeros(1, length(ifor));
        
        get_segtens           
        rdiff = get_rdiff(rbead,1,ipt);
        temp = circshift(rdiff,1,2);
        segtens_scalar = -sum(temp.*segtens);

        beadidx_b = zeros(1,length(segtens));
        for jseg = 1:length(segtens)
            if any(ifor==jseg)
                beadidx_b(jseg) = 1;
            else
                beadidx_b(jseg) = beadidx_b(jseg-1) + 1;
            end
        end
        
        for ifil = 1:length(ifor)
            ix1 = find(any(ifor(ifil):ipt(ifil) == find(sum(xmat))'), 2);
            if(ix1 <= 4)
                isx(ifil) = 1;
            end
            inode = find(fmmat(ifil, :));
            if sum(fmmat(:, inode)) ~= 1
                segtens_scalar(ifor(ifil):ipt(ifil)) = nan;
            end
        end

        for i = 2:bdmax
            tseg(i-1) = nanmean(segtens_scalar(beadidx_b==i));
            tseg_sd(i-1) = nanstd(segtens_scalar(beadidx_b==i));
        end

        tseg_store(idx, :) = tseg;
        
        phi = cart2pol(rbead(1, ifor), rbead(2, ifor));
        phi_hat = [-sin(phi); cos(phi); 0*phi];
        pol = sign(sum(phi_hat.*rdiff(:, ifor)));
        vf_store = [vf_store, sum(phi_hat.*vbead(:, ifor)).*pol];
        tb_store = [tb_store segtens_scalar(beadidx_b==2)];
        l_store = [l_store (ipt-ifor)*0.1];
        isx_store = [isx_store isx];
    end

    if(isempty(vf_store))
        disp('No data with these parameters')
        continue
    end
    
    %% tension v. filament length
    clear t_bin tsd_bin tsem_bin
    lbin = 0:0.1:10;
    for ib = 1:(length(lbin)-1)
        t_bin(ib) = nanmean(tb_store(lbin(ib) <= l_store & l_store < lbin(ib+1)));
        tsd_bin(ib) = nanstd(tb_store(lbin(ib) <= l_store & l_store < lbin(ib+1)));
        tsem_bin(ib) = tsd_bin(ib)/sqrt(length(tb_store(lbin(ib) <= l_store & l_store < lbin(ib+1)))-1);
    end
    % figure
%     errorbar(x, nanmean(tseg_store), nanstd(tseg_store), 'displayname', ['\gamma = ' num2str(gmlist(i0)) 'pN s/\mum'])
%     tb_store = [tb_store; tseg_store(:, 1)];
%     figure
%     histogram(tb_store, -20:2:60)
    figure
    hold on
    plot(l_store, tb_store, 'r.')
    errorbar(lbin(1:end-1), t_bin, tsem_bin, 'ko', 'MarkerFaceColor', 'auto')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    ylabel('Barbed end tension (pN)')
    xlabel('Filament length (\mum)')
    title(['i0=' num2str(i0)])
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
%     xlim([0 5])
    ylim([quantile(tb_store, 0.01), quantile(tb_store, 0.99)])
    saveas(gcf, [pfx '_ltbarb_i' num2str(i0) '.png']);
    saveas(gcf, [pfx '_ltbarb_i' num2str(i0) '.fig']);

    %% tension v. velocity
    clear t_bin tsd_bin tsem_bin
    vfbin = -0.2:0.005:0.2;
    for ib = 1:(length(vfbin)-1)
        t_bin(ib) = nanmean(tb_store(vfbin(ib) <= vf_store & vf_store < vfbin(ib+1)));
        tsd_bin(ib) = nanstd(tb_store(vfbin(ib) <= vf_store & vf_store < vfbin(ib+1)));
        tsem_bin(ib) = tsd_bin(ib)/sqrt(length(tb_store(vfbin(ib) <= vf_store & vf_store < vfbin(ib+1)))-1);
    end
    figure
    hold on
    plot(vf_store*1e3, tb_store, 'r.')
    errorbar(vfbin(1:end-1)*1e3, t_bin, tsem_bin, 'ko', 'MarkerFaceColor', 'auto')
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
%     ylim([-20 60])
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    xlabel('Node velocity (nm/s)')
    ylabel('Barbed end tension (pN)')
    title(['i0=' num2str(i0)])
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
    xlim([quantile(vf_store, 0.01), quantile(vf_store, 0.99)]*1e3)
    ylim([quantile(tb_store, 0.01), quantile(tb_store, 0.99)])
    saveas(gcf, [pfx '_vtbarb_i' num2str(i0) '.png']);
    saveas(gcf, [pfx '_vtbarb_i' num2str(i0) '.fig']);

%     figure
%     hold on
%     plot(l_store, vf_store*1e3, 'ko')
%     plot(xlim, [0 0], 'k--', 'linewidth', 2)
%     set(gca, 'Linewidth', 3)
%     set(gca, 'FontSize', 24)
%     set(gca, 'Position', [0.2 0.2 0.7 0.7]);
%     xlabel('Filament length (\mum)')
%     ylabel('Node velocity (nm/s)')
%     title(['i0=' num2str(i0)])
%     plot(xlim, [0 0], 'k--', 'linewidth', 2)
%     saveas(gcf, [pfx '_vlfil_i' num2str(i0) '.png']);
%     saveas(gcf, [pfx '_vlfil_i' num2str(i0) '.fig']);

    %% velocity v. length
    lbin = 0:0.1:10;
    for ib = 1:(length(lbin)-1)
        vm_bin(ib) = 1e3*nanmean(vf_store(lbin(ib) <= l_store & l_store < lbin(ib+1)));
        vsd_bin(ib) = 1e3*nanstd(vf_store(lbin(ib) <= l_store & l_store < lbin(ib+1)));
        vsem_bin(ib) = vsd_bin(ib)/length(vf_store(lbin(ib) <= l_store & l_store < lbin(ib+1))-1);
    end
    figure
    hold on
    plot(l_store, vf_store*1e3, 'r.')
    errorbar(lbin(1:end-1), vm_bin, vsem_bin, 'ko', 'MarkerFaceColor', 'auto')
%     vlfit = fit(lbin(~isnan(vm_bin))', vm_bin(~isnan(vm_bin))', 'poly1');
%     vlfit = fit(l_store', vf_store', 'poly1');
%     plot(vlfit)
%     vl_slope(irun) = vlfit.p1;
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    xlabel('Filament length (\mum)')
    ylabel('Node velocity (nm/s)')
    title(['i0=' num2str(i0)])
    ylim([quantile(vf_store, 0.01), quantile(vf_store, 0.99)]*1e3)
    saveas(gcf, [pfx '_lvf_i' num2str(i0) '.png']);
    saveas(gcf, [pfx '_lvf_i' num2str(i0) '.fig']);
    
    pcomp(i0) = sum(tb_store < 0)/sum(~isnan(tb_store));
    pback(i0) = sum(vf_store < 0)/sum(~isnan(vf_store));
end
