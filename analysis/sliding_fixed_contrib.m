clear 
tfixed = zeros(5);
tcirc = zeros(5);
tslide = zeros(5);
tlap = zeros(5);
pval = zeros(5);
pave = zeros(5);
vpagg = [];
vmagg = [];
vforagg = [];
vagg = [];
dt = 0.2;
nstep = 1;

% for irun = 4:5:25 % wt runs
for irun = 1:25 % all runs
    clearvars -except irun dt nstep ttotal tcirc tlap tfixed tslide ...
        pval pave vpagg vmagg vforagg vagg 
    load(['gm_2_' num2str(irun) '_00sec.mat']);
    load(['gm_2_' num2str(irun) '_600sec.mat']);
    rk1_v32
    tension_laplace_v1
    tlap(irun) = tens;
    tension_circ_v4
    tcirc(irun) = tens;
%     bmmat = bm_nbr(rbead, rmyo, ifor, ipt, rcapmyo_long, rcapmyo_short, fmmat);
%     bpmat = bp_nbr(rbead, rmyp, ifor, ipt, rcapmyp);
    %% polarity distribution
    p = zeros(1, length(rmyo));
    v = zeros(1, length(rmyo));
    for inode = 1:length(rmyo) % loop over all nodes
        p(inode) = 0;
        % local phihat direction, circumferential around ring
        phi = atan2(rmyo(2, inode), rmyo(1, inode));
        phihat = [-sin(phi), cos(phi), 0];
        v(inode) = phihat*vmyo(:, inode);

        %% calculating polarity by total filament length
        if(sum(fmmat(:, inode)) == 1)
        for idx_for = ifor(find(fmmat(:, inode))) % loop over fils on this node
            % find associated pointed end
            idx_pt = ipt(find(ifor == idx_for));
            % tangent vector
            t = rbead(:, idx_pt) - rbead(:, idx_for);
            t = t/sqrt(sum(t.^2));
            % add/subtract length of filament depending on polarity
            p(inode) = p(inode) + sign(phihat*t)*(idx_pt - idx_for);
        end
        end
        % take sign to determine winning polarity
        p(inode) = sign(p(inode));
        %% calculating polarity by velocity
%         p(inode) = sign(v(inode));
    end
    
    vpagg = [vpagg v(p==1)];
    vmagg = [vmagg v(p==-1)];
    vforagg = [vforagg v(p~=0)];
    vagg = [vagg v];
    
%     figure
%     hold on
%     bins = -0.15:0.002:0.15;
%     histogram(v, bins);
%     histogram(v(~(p ==0)), bins);
% %     histogram(v(p == -1), bins)
% %     histogram(v(p == 1), bins)
%     pval(irun) = ranksum(v(p == 1), v(p==-1));
    
%     figure
%     histogram(p);
%     pave(irun) = mean(p);
%     [junk pval(irun)] = ttest(p);
    [junk, pval(irun)] = ttest2(v(p == 1), v(p == -1));
    
    fnode = zeros(length(rmyo));
    % stall force is renormalized due to cycling
    sfpf = stall_force_pf(bmmat, fone, nsat, nhead, fhead);
    % stall force is renormalized due to cycling
    sfpf_myp = stall_force_pf(bpmat, fone, nsat, nhead, fheadmyp);
    %% tension contribution
    for ibead = 1:length(rbead) % loop over all actin beads
        ifil = find(ipt >= ibead, 1); % identify host filament
        inode = find(fmmat(ifil, :)); % identify host node
        node_list = find(bmmat( ibead, :)); % find nodes pulling this bead
        x = (ibead - ifor(ifil))*dbead; % distance from barbed end
        lfil = (ipt(ifil) - ifor(ifil))*dbead;
        
        for jnode = node_list % loop over guest nodes
            % tangent vector
            t = rbead(:, ibead) - rbead(:, ibead + 1);
            t = t/sqrt(sum(t.^2));
            if(p(jnode) == p(inode))
                v = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
                tfixed(irun) = tfixed(irun) + x*sfpf(jnode) * ...
                    (1 - t'*v/vmyo0);
            else
                v = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
                tslide(irun) = tslide(irun) + x*sfpf(jnode) * ...
                    (1 - t'*v/vmyo0);
            end
        end
        
        myp_list = find(bpmat( ibead, :));
        for imyp = myp_list
            % !!!WARNING: MYP2 VELOCITY NOT SAVED; ASSUMED ZERO!!!
            v = vmyp(:, imyp) - vbead(:, ibead);
            if(p(inode) == 0)
                tfixed(irun) = tfixed(irun) + x*sfpf_myp(imyp) * ...
                    (1 - t'*v/vmyo0);
            else
                tslide(irun) = tslide(irun) + x*sfpf_myp(imyp) * ...
                    (1 - t'*v/vmyo0);
            end
        end
    end
    tfixed(irun) = tfixed(irun)/(2*pi*r_ring);
    tslide(irun) = tslide(irun)/(2*pi*r_ring);
end

gm = [0.2 0.5 0.8 1 1.5];
figure
hold on
errorbar(gm, mean(tfixed'), std(tfixed'), 'ro', 'linewidth', 2)
errorbar(gm, mean(tslide'), std(tslide'), 'bo', 'linewidth', 2)
errorbar(gm, mean(tfixed') + mean(tslide'), std(tfixed') + std(tslide'), ...
    'ko', 'linewidth', 2)
errorbar(gm, mean(tlap'), std(tlap'), 'mo', 'linewidth', 2)
errorbar(gm, mean(tcirc'), std(tcirc'), 'co', 'linewidth', 2)
legend('fixed', 'sliding', 'total', 'laplace')
ylabel('Tension (pN)')
xlabel('Relative node anchor drag coefficient')
set(gca, 'Linewidth', 3)
set(gca, 'FontSize', 24)
set(gca, 'Position', [0.2 0.2 0.7 0.7]);