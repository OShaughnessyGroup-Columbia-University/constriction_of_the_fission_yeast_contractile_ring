tfixed = zeros(5);
tslide = zeros(5);
pval = zeros(5);
pave = zeros(5);
vpagg = [];
vmagg = [];
vforagg = [];
vagg = [];

load(['gm_2_1_00sec.mat']);
for irun = 1:25
    load(['gm_2_' num2str(irun) '_600sec.mat']);
    bmmat = bm_nbr(rbead, rmyo, ifor, ipt, rcapmyo_long, rcapmyo_short, fmmat);
    bpmat = bp_nbr(rbead, rmyp, ifor, ipt, rcapmyp);
    %% polarity distribution
    p = zeros(1, length(rmyo));
    v = zeros(1, length(rmyo));
    for inode = 1:length(rmyo) % loop over all nodes
        p(inode) = 0;
        % local phihat direction, circumferential around ring
        phi = atan2(rmyo(2, inode), rmyo(1, inode));
        phihat = [-sin(phi), cos(phi), 0];
        v(inode) = phihat*vmyo(:, inode);
        
        for idx_for = ifor(find(fmmat(:, inode))) % loop over fils on this node
            % tangent vector
            t = rbead(:, idx_for) - rbead(:, idx_for + 1);
            t = t/sqrt(sum(t.^2));
            % find associated pointed end
            idx_pt = ipt(find(ifor == idx_for));
            % add/subtract length of filament depending on polarity
            p(inode) = p(inode) + sign(phihat*t)*(idx_pt - idx_for);
        end
        % take sign to determine winning polarity
        p(inode) = sign(p(inode));
%         if p(inode) < -5
%             p(inode) = -1;
%         elseif -5 <= p(inode) && p(inode) < 5
%             p(inode) = 0;
%         else
%             p(inode) = 1;
%         end
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
    [junk pval(irun)] = ttest2(v(p == 1), v(p == -1));
    
    fnode = zeros(length(rmyo));
    %% tension contribution
    for ifil = 1:length(ifor) % loop over filaments
        lfil = (ipt(ifil) - ifor(ifil))*0.1; % length of filament
        inode = find(fmmat(ifil, :)); % identify host node
        
        for ibead = ifor(ifil):ipt(ifil) % loop over beads in this fil
            node_list = find(bmmat( ibead, :)); % find nodes pulling this bead
            for jnode = node_list % loop over guest nodes
                % tangent vector
                t = rbead(:, ibead) - rbead(:, ibead + 1);
                t = t/sqrt(sum(t.^2));
                if(p(jnode) == p(inode))
                    v = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
                    tfixed(irun) = tfixed(irun) + fhead * (1 - t'*v/vmyo0);
                else
                    v = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
                    tslide(irun) = tslide(irun) + fhead * (1 - t'*v/vmyo0);
                end
            end  
            
            myp_list = find(bpmat( ibead, :));
            for imyp = myp_list
                % !!!WARNING: DON'T HAVE MYP2 VELOCITY; ASSUMED ZERO!!!
                v = zeros(3, 1) - vbead(:, ibead);
                if(p(inode) == 0)
                    tfixed(irun) = tslide(irun) + fheadmyp * (1 - t'*v/vmyo0);
                else
                    tslide(irun) = tslide(irun) + fheadmyp * (1 - t'*v/vmyo0);
                end
            end
        end
        
    end
    tfixed(irun) = tfixed(irun)*length(rbead)*0.1/(2*pi*r_ring) / length(ifor);
    tslide(irun) = tslide(irun)*length(rbead)*0.1/(2*pi*r_ring) / length(ifor);
end

gm = [0.2 0.5 0.8 1 1.5];
figure
hold on
errorbar(gm, mean(tfixed'), std(tfixed'), 'ro')
errorbar(gm, mean(tslide'), std(tslide'), 'bo')
errorbar(gm, mean(tfixed') + mean(tslide'), std(tfixed') + std(tslide'), 'ko')