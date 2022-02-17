tfixed = 0;
tslide = 0;
vave = 0;

% for irun = 4:5:25 % wt runs
bmmat = bm_nbr(rbead, rmyo, ifor, ipt, rcapmyo_long, rcapmyo_short, fmmat);
bpmat = bp_nbr(rbead, rmyp, ifor, ipt, rcapmyp);
%% polarity distribution
ptot = zeros(1, length(rmyo));
p = zeros(2, length(rmyo));
v = zeros(1, length(rmyo));
for inode = 1:length(rmyo) % loop over all nodes
    ptot(inode) = 0;
    % local phihat direction, circumferential around ring
    phi = atan2(rmyo(2, inode), rmyo(1, inode));
    phihat = [-sin(phi), cos(phi), 0];
    v(inode) = phihat*vmyo(:, inode);
    
    %% calculating polarity by total filament length
    for idx_for = ifor(find(fmmat(:, inode))) % loop over fils on this node
        if(isgood(find(ifor == idx_for), rbead, rmyo, ifor, ipt, r_ring))
            % find associated pointed end
            idx_pt = ipt(find(ifor == idx_for));
            % tangent vector
            t = rbead(:, idx_for+1) - rbead(:, idx_for);
            t = t/sqrt(sum(t.^2));
            % increment +/- filament count
            if(sign(phihat*t) == 1)
                p(1, inode) = p(1, inode) + 1;
            elseif(sign(phihat*t) == -1)
                p(2, inode) = p(2, inode) + 1;
            else
                disp('null filament!')
            end
            % add/subtract length of filament depending on polarity
            ptot(inode) = ptot(inode) + sign(phihat*t)*(idx_pt - idx_for);
%         figure
%         disp(['p = ' num2str(sign(phihat*t))])
%         plot3(rbead(1, idx_for), ...
%             rbead(2, idx_for),...
%             rbead(3, idx_for), 'o')
%         hold on
%         plot3(rbead(1, idx_for:idx_pt), ...
%             rbead(2, idx_for:idx_pt), ...
%             rbead(3, idx_for:idx_pt))
%         view(gca, [0, 0, 100])
%         xlim([-2 2])
%         ylim([-2 2])
%         close
        end
    end
    % take sign to determine winning polarity
    ptot(inode) = sign(ptot(inode));
%     %% calculating polarity by velocity
%     p(inode) = sign(v(inode));
end
vave = mean(abs(v(p(1, :)~=0 | p(2, :)~=0)));

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
        if(ptot(jnode) == ptot(inode))
            dv = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
            tfixed = tfixed + x*sfpf(jnode) * ...
                (1 - t'*dv/vmyo0);
        else
            dv = vmyo(:, jnode) - vbead(:, ibead); % relative velocity
            tslide = tslide + x*sfpf(jnode) * ...
                (1 - t'*dv/vmyo0);
        end
    end
    
    myp_list = find(bpmat( ibead, :));
    for imyp = myp_list
        if(~isempty(vmyp))
            dv = vmyp(:, imyp) - vbead(:, ibead);
        else
            dv = - vbead(:, ibead);
        end
        
        if(ptot(inode) == 0)
            tfixed = tfixed + x*sfpf_myp(imyp) * ...
                (1 - t'*dv/vmyo0);
        else
            tslide = tslide + x*sfpf_myp(imyp) * ...
                (1 - t'*dv/vmyo0);
        end
    end
end
tfixed = tfixed/(2*pi*r_ring) ; %% SOME DOUBLE COUNTING; UNEXPLAINED
tslide = tslide/(2*pi*r_ring) ; %% SOME DOUBLE COUNTING; UNEXPLAINED