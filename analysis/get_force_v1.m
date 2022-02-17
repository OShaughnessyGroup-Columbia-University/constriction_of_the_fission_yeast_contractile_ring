% change from original: use rexc_n for node-node excluded volume

function [force, gbeadfv, gppf, bmmatpt, bpmatpt] = get_force_v1(rbead,rmyo,rmyp,rdiff,ifor,ipt,bancm,bmmat,bpmat,nhead,fhead,fheadmyp,fone,nsat,rcapmyo_u, rcapmyo_long, rcapmyo_short,rcapmyp,kcap,kcap_p,dbead,vmyo0,kod2,r_ring,kwall,kexc,rexc,rexc_n,kex_act,rex_act,xmat,rx0,kx,d_cdc15,d_myo,geometry)
    % calculate the forces on all particles given the velocities
    
    % will be repeatedly run.
    
    % myo2 pulling stall force
        % stall force per filament of each myosin cluster, full vector
        sfpf_myo = stall_force_pf (bmmat,fone,nsat,nhead,fhead);
        % total stall force amplitude on each bead, full vector
        fsamp = bmmat * sfpf_myo;
        % make it a row vector
        fsamp = fsamp';
        % the contribution from myosin pulling force
        if isempty(fsamp)
            fbead = zeros(size(rbead));
        else
            fbead = rdiff .* [fsamp; fsamp; fsamp];
        end
        % sum of rdiff of all beads that interact with each myosin cluster
        rdiffsum = rdiff * bmmat;
            rdiffp = rdiff;
            rdiffp(:,ipt) = rdiffp(:,ipt-1);            % ridffp is the unit tangent vector at all beads including pointed ends
        % the contribution from myosin pulling force
        if isempty(sfpf_myo)
            fmyo = double.empty(3,0);
        elseif isempty(rdiffsum)
            fmyo = zeros(size(rmyo));
        else
            fmyo = - rdiffsum .* repmat(sfpf_myo',3,1);
        end

        % effective drag coefficient (gamma_pull)of each myosin cluster. NOTE: PER FILAMENT!!!
        gppf = sfpf_myo / vmyo0;
        % the sum of effective drags from all myosin clusters
        if ~isempty(rbead)&&isempty(bmmat)
            gbeadfv = zeros(size(rbead,2),1);
        else
            gbeadfv = bmmat * gppf;
        end
        % make it a row vector
        gppf = gppf';

    % myp2 pulling stall force
        % stall force per filament of each myosin cluster, full vector
        sfpf_myp = stall_force_pf (bpmat,fone,nsat,nhead,fheadmyp);
        % total stall force amplitude on each bead, full vector
        fsamp = bpmat * sfpf_myp;
        % make it a row vector
        fsamp = fsamp';
        % the contribution from myosin pulling force
        if ~isempty(fsamp)
            fbead = fbead + rdiff .* [fsamp; fsamp; fsamp];
        end
        % sum of rdiff of all beads that interact with each myosin cluster
        rdiffsum = rdiff * bpmat;
            rdiffp = rdiff;
            rdiffp(:,ipt) = rdiffp(:,ipt-1);            % ridffp is the unit tangent vector at all beads including pointed ends
        % the contribution from myosin pulling force
        if isempty(sfpf_myp)
            fmyp = double.empty(3,0);
        elseif isempty(rdiffsum)
            fmyp = zeros(size(rmyp));
        else
            fmyp = - rdiffsum .* repmat(sfpf_myp',3,1);
        end

        % effective drag coefficient (gamma_pull)of each myosin cluster. NOTE: PER FILAMENT!!!
        gppfp = sfpf_myp / vmyo0;
        % the sum of effective drags from all myosin clusters
        if ~isempty(bpmat)
            gbeadfv = gbeadfv + bpmat * gppfp;
        end
        % make it a row vector
        gppfp = gppfp';
        gppf = [gppf, gppfp];

    % myo2 capture force  
        % use bmmatpt which includes pointed end
        bmmatpt = ptcap(rbead, rmyo, ipt, bancm, rcapmyo_u, rcapmyo_long, rcapmyo_short, bmmat);
        [ib, im] = find(bmmatpt);
        for i = 1:numel(ib)
            this_ib = ib(i);
            this_im = im(i);
            rbm = rmyo(:,this_im) - rbead(:,this_ib);   % r from b to m
            rbm_perp = rbm - sum(rbm .* rdiffp(:,this_ib)) * rdiffp(:,this_ib);
            f_m2b = kcap * rbm_perp;
            fbead(:,this_ib) = fbead(:,this_ib) + f_m2b;
            fmyo(:,this_im) = fmyo(:,this_im) - f_m2b;
        end
        
    % myp2 capture force  
        % use bpmatpt which includes pointed end
        bpmatpt = ptcap_myp(rbead, rmyp, ipt, rcapmyp, bpmat);
        [ib, im] = find(bpmatpt);
        for i = 1:numel(ib)
            this_ib = ib(i);
            this_im = im(i);
            rbm = rmyp(:,this_im) - rbead(:,this_ib);   % r from b to m
            rbm_perp = rbm - sum(rbm .* rdiffp(:,this_ib)) * rdiffp(:,this_ib);
            f_m2b = kcap_p * rbm_perp;
            fbead(:,this_ib) = fbead(:,this_ib) + f_m2b;
            fmyp(:,this_im) = fmyp(:,this_im) - f_m2b;
        end
    
        %%
    % actin filament bending force
        % construct fi
        fi = true(1,size(rbead,2));
        fi(ipt) = false;
        fi(ifor) = false;
        fi = find(fi);
        % construct f_i-1
        fimo = true(1,size(rbead,2));
        fimo(ifor) = false;
        fimo(ifor + 1) = false;
        fimo = find(fimo);
        % construct f_i+1
        fipo = true(1,size(rbead,2));
        fipo(ipt) = false;
        fipo(ipt - 1) = false;
        fipo = find(fipo);
        % construct t_i-1, t_i+1 and t_i+2
        tim1 = circshift(rdiff,[0 2]);
        ti = circshift(rdiff,[0 1]);
        tip1 = rdiff;
        tip2 = circshift(rdiff,[0 -1]);
        % bending force
        fb = zeros(size(fbead));
        for i = fimo
            fb(:,i) = kod2 * (tim1(:,i) - (tim1(:,i)' * ti(:,i)) * ti(:,i));
        end
        for i = fi
            fb(:,i) = fb(:,i) + kod2 * (-(1 + ti(:,i)' * tip1(:,i)) * ti(:,i)...
                + (ti(:,i)' * tip1(:,i) + 1) * tip1(:,i));
            if any(ifor+1==i)   % the first segment after a formin has a different length
                temp = 100*dbead;%norm(rbead(:,i)-rbead(:,i-1));       % length of this segment
                fb(:,i) = fb(:,i) - kod2 * (-ti(:,i)' * tip1(:,i)) * ti(:,i) + kod2 * (-ti(:,i)' * tip1(:,i)) * ti(:,i) * dbead / temp;
                fb(:,i) = fb(:,i) - kod2 * tip1(:,i) + kod2 * tip1(:,i) * dbead / temp;
            end
        end
        for i = fipo
            fb(:,i) = fb(:,i) + kod2 * ((tip1(:,i)' * tip2(:,i)) * tip1(:,i) - tip2(:,i));
            if any(ifor+1==i)   % the first segment after a formin has a different length
                temp = 100*dbead;%norm(rbead(:,i)-rbead(:,i-1));       % length of this segment
                fb(:,i) = fb(:,i) - kod2 * ((tip1(:,i)' * tip2(:,i)) * tip1(:,i) - tip2(:,i))...
                    + kod2 * ((tip1(:,i)' * tip2(:,i)) * tip1(:,i) - tip2(:,i)) * dbead / temp;
            end
        end
% % %         % a maximum of maxfb is applied
% % %         fbnorm = sqrt(sum(fb .* fb));
% % %         ind = fbnorm > maxfb;
% % %         for i = 1:3
% % %             fb(i,ind) = fb(i,ind) * maxfb ./ fbnorm(ind);
% % %         end
        fbead = fbead + fb;
    %% plasma membrane -- restricting everything inside
    if geometry == 'sph'
        % decide which beads / myo go out of the membrane
        radbead = sqrt(sum(rbead .* rbead));
        rring = r_ring;
        %rring = sqrt(rring^2 + .05^2)
        outind = find(radbead > rring);
        for i = outind
            fwall = - kwall * rbead(:,i) / radbead(i) * (radbead(i) - rring);
            fbead(:,i) = fbead(:,i) + fwall;
        end

        radmyo = sqrt(sum(rmyo .* rmyo));
        outind = find(radmyo > rring);
        for i = outind
            fwall = - kwall * rmyo(:,i) / radmyo(i) * (radmyo(i) - rring);
            fmyo(:,i) = fmyo(:,i) + fwall;
        end

        radmyp = sqrt(sum(rmyp .* rmyp));
        outind = find(radmyp > rring);
        for i = outind
            fwall = - kwall * rmyp(:,i) / radmyp(i) * (radmyp(i) - rring);
            fmyp(:,i) = fmyp(:,i) + fwall;
        end
    end
    
    if strcmp(geometry,'cyl')&&~isempty(rbead)&&~isempty(rmyo)
        % decide which beads / myo go out of the membrane
        radbead = sqrt(sum(rbead(1:2,:) .* rbead(1:2,:)));
        rring = r_ring;
        %rring = sqrt(rring^2 + .05^2)
        outind = find(radbead > rring);
        for i = outind
            fwall = - kwall * rbead(1:2,i) / radbead(i) * (radbead(i) - rring);
            fbead(1:2,i) = fbead(1:2,i) + fwall;
        end

        radmyo = sqrt(sum(rmyo(1:2,:) .* rmyo(1:2,:)));
        outind = find(radmyo > rring);
        for i = outind
            fwall = - kwall * rmyo(1:2,i) / radmyo(i) * (radmyo(i) - rring);
            fmyo(1:2,i) = fmyo(1:2,i) + fwall;
        end

        radmyp = sqrt(sum(rmyp(1:2,:) .* rmyp(1:2,:)));
        outind = find(radmyp > rring-rcapmyp);
        for i = outind
            fwall = - kwall * rmyp(1:2,i) / radmyp(i) * (radmyp(i) - rring+rcapmyp);
            fmyp(1:2,i) = fmyp(1:2,i) + fwall;
        end
    end
    %% myosin excluded volume
        % myo2-myo2
        for i = 1:size(rmyo,2)
            % difference between the x component of all myosin clusters and
            % the x component of this myosin cluster
            delx = rmyo(1,:) - rmyo(1,i);
            % these are the ones with |delx| < rexc
            ind = abs(delx) < rexc_n;
            % exclude itself
            ind(i) = false;
            % switch from boolean mode to index mode
            ind = find(ind);
            % difference in y direction 
            dely = rmyo(2,ind) - rmyo(2,i);
            % these are the ones with |delx| < rexc and |dely| < rexc
            ind = ind(abs(dely) < rexc_n);
            % if no myosin satisfies this criterion, skip to next
            if isempty(ind)
                continue
            end
            % difference in all 3 directions
            delx = delx(ind);
            dely = dely(abs(dely) < rexc_n);
            delz = rmyo(3,ind) - rmyo(3,i);
            % total distance 
            delr = delx .* delx + dely .* dely + delz .* delz;
            delr = sqrt(delr);
            % these are the ones with delr < rexc
            ind = ind(delr < rexc_n);
            % if no other myosins within rexc, skip to next 
            if isempty(ind)
                continue
            end
            % delx, dely, delz and delr for those
            delx = delx(delr < rexc_n);   
            dely = dely(delr < rexc_n);
            delz = delz(delr < rexc_n);
            delr = delr(delr < rexc_n);
            
            % go through every myosin within rexc
            for j = 1:length(ind)
                % magnitude of force
                fmag = kexc * (rexc_n - delr(j));
                % unit vector of delr
                delrunit = [delx(j); dely(j); delz(j)];
                delrunit = delrunit / norm(delrunit);
                % vector force
                f = - fmag * delrunit;
                % incoporate into fmyo
                fmyo(:,i) = fmyo(:,i) + f;
            end
        end   
        
        % myp2-myp2
        for i = 1:size(rmyp,2)
            % difference between the x component of all myosin clusters and
            % the x component of this myosin cluster
            delx = rmyp(1,:) - rmyp(1,i);
            % these are the ones with |delx| < rexc
            ind = abs(delx) < rexc;
            % exclude itself
            ind(i) = false;
            % switch from boolean mode to index mode
            ind = find(ind);
            % difference in y direction 
            dely = rmyp(2,ind) - rmyp(2,i);
            % these are the ones with |delx| < rexc and |dely| < rexc
            ind = ind(abs(dely) < rexc);
            % if no myosin satisfies this criterion, skip to next
            if isempty(ind)
                continue
            end
            % difference in all 3 directions
            delx = delx(ind);
            dely = dely(abs(dely) < rexc);
            delz = rmyp(3,ind) - rmyp(3,i);
            % total distance 
            delr = delx .* delx + dely .* dely + delz .* delz;
            delr = sqrt(delr);
            % these are the ones with delr < rexc
            ind = ind(delr < rexc);
            % if no other myosins within rexc, skip to next 
            if isempty(ind)
                continue
            end
            % delx, dely, delz and delr for those
            delx = delx(delr < rexc);   
            dely = dely(delr < rexc);
            delz = delz(delr < rexc);
            delr = delr(delr < rexc);
            
            % go through every myosin within rexc
            for j = 1:length(ind)
                % magnitude of force
                fmag = kexc * (rexc - delr(j));
                % unit vector of delr
                delrunit = [delx(j); dely(j); delz(j)];
                delrunit = delrunit / norm(delrunit);
                % vector force
                f = - fmag * delrunit;
                % incoporate into fmyo
                fmyp(:,i) = fmyp(:,i) + f;
            end
        end
        % myo2-myp2
        for i = 1:size(rmyp,2)
            % difference between the x component of all myo2 clusters and
            % the x component of this myp2 cluster
            delx = rmyo(1,:) - rmyp(1,i);
            % these are the ones with |delx| < rexc
            ind = abs(delx) < rexc;
            % switch from boolean mode to index mode
            ind = find(ind);
            % difference in y direction 
            dely = rmyo(2,ind) - rmyp(2,i);
            % these are the ones with |delx| < rexc and |dely| < rexc
            ind = ind(abs(dely) < rexc);
            % if no myosin satisfies this criterion, skip to next
            if isempty(ind)
                continue
            end
            % difference in all 3 directions
            delx = delx(ind);
            dely = dely(abs(dely) < rexc);
            delz = rmyo(3,ind) - rmyp(3,i);
            % total distance 
            delr = delx .* delx + dely .* dely + delz .* delz;
            delr = sqrt(delr);
            % these are the ones with delr < rexc
            ind = ind(delr < rexc);
            % if no other myosins within rexc, skip to next 
            if isempty(ind)
                continue
            end
            % delx, dely, delz and delr for those
            delx = delx(delr < rexc);   
            dely = dely(delr < rexc);
            delz = delz(delr < rexc);
            delr = delr(delr < rexc);
            
            % go through every myo2 within rexc
            for j = 1:length(ind)
                % magnitude of force
                fmag = kexc * (rexc - delr(j));
                % unit vector of delr
                delrunit = [delx(j); dely(j); delz(j)];
                delrunit = delrunit / norm(delrunit);
                % vector force
                f = - fmag * delrunit;
                % incoporate into fmyo
                fmyp(:,i) = fmyp(:,i) + f;
                fmyo(:,ind(j)) = fmyo(:,ind(j)) - f;
            end
        end
        %% node-myp2 excluded volume
    % go through every myo2 cluster
        for i = 1:size(rmyo,2)
            % find location of cdc15
            this_rmyo = rmyo(:,i);
            r_cdc15 = this_rmyo * (r_ring-d_cdc15)/(r_ring-d_myo);
            r_cdc15(3) = this_rmyo(3);
            % difference between the x component of all myp2 clusters and
            % the x component of this cdc15
            delx = rmyp(1,:) - r_cdc15(1);
            % these are the ones with |delx| < rexc/2. Note: rexc is the
            % excluded diameter! Unfortunate name.
            ind = abs(delx) < rexc;
            % exclude itself
            ind(i) = false;
            % switch from boolean mode to index mode
            ind = find(ind);
            % difference in y direction 
            dely = rmyp(2,ind) - r_cdc15(2);
            % these are the ones with |delx| < rexc and |dely| < rexc
            ind = ind(abs(dely) < rexc);
            % if no myosin satisfies this criterion, skip to next
            if isempty(ind)
                continue
            end
            % difference in all 3 directions
            delx = delx(ind);
            dely = dely(abs(dely) < rexc);
            delz = rmyp(3,ind) - r_cdc15(3);
            % total distance 
            delr = delx .* delx + dely .* dely + delz .* delz;
            delr = sqrt(delr);
            % these are the ones with delr < rexc
            ind = ind(delr < rexc);
            % if no other myosins within rexc, skip to next 
            if isempty(ind)
                continue
            end
            % delx, dely, delz and delr for those
            delx = delx(delr < rexc);   
            dely = dely(delr < rexc);
            delz = delz(delr < rexc);
            delr = delr(delr < rexc);
            
            % go through every myp2 within rexc/2
            for j = 1:length(ind)
                % magnitude of force
                fmag = kexc * (rexc - delr(j));
                % unit vector of delr
                delrunit = [delx(j); dely(j); delz(j)];
                delrunit = delrunit / norm(delrunit);
                % vector force
                f = - fmag * delrunit;
                % incoporate into fmyo
                fmyo(:,i) = fmyo(:,i) + f;
                % incorporate into fmyp
                fmyp(:,ind(j)) = fmyp(:,ind(j)) - f;
            end
        end 
        %% actin bead-based excluded volume
        temp = act_excl_force(rbead,kex_act,rex_act,2*rex_act,ifor);
        fbead = fbead + temp;
        %% crosslinker elastic force
        
        [row, col] = find(xmat);
        for i = 1:numel(row)
            thisrow = row(i);
            thiscol = col(i);
            drv = rbead(:,thisrow) - rbead(:,thiscol);      % vector of dr
            dra = sqrt(sum(drv .* drv));                    % amplitude of dr
            fspa = kx *(rx0 - dra);                         % amplitude of spring force. Positive means compression.
            fspv = fspa * drv / dra;                  % vector of spring force on the bead thisrow.
            fbead(:,thisrow) = fbead(:,thisrow) + fspv;
            fbead(:,thiscol) = fbead(:,thiscol) - fspv;         
        end
        
%%
        % concatenate fbead and fmyo into a 1 by 3*(nbead + nmyo) vector
        force = [fbead, fmyo, fmyp];
        force = force(:)';
end