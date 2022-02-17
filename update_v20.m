% change from v4: fixed the bug that causes xlinker number to decrease in time
% change from v5: use rangesearch in xlinker finding to speed up.
% change from v6: do not remove Myp2p that is not bound to any actin
% (for parameter scanning purpose)
% change from v7: formin nucleates actin in such a way that it avoids
% overlapping with existing actin
% v8 skipped
% change from v9: take dt as input, no myosin rebinding (to be called at
% every time step), corrected an error in calculating the distance from proposed new
% actin bead to existing beads
% change from v10: also use new formin position if can't find a place to
% insert actin; a bug corrected such that a new actin bead does not need to
% avoid existing formin other than its own formin; a bug corrected
% regarding beads polymerizing towards membrane; 
% change from v11: if an actin segment becomes too long, remove the
% filament
% change from v12: myo2 avoids existing myo2 when binding
% change from v13: get rid of myo2 that fly too far away
% change from v14: myp2 binds with homogeneous spacial probability to
% actin, restore unbinding of free myp2; actin segment abnormal length
% threshold changed to 2*dbead
% change from v15: maintain formin number

function [rbead,rmyo,rmyp,xmat,fmmat,ipt,ifor,bancf,bancm,dbead_first,bpmat_trial] = update_v20(rbead,rmyo,rmyp,bpmat,xmat,rhoain,rxbind,rex_act,rexc_n,fmmat,bancf,bancm,ifor,ipt,kofffor,...
    koffmyo,koffmyp,dt,d_for,d_myo,dbead,dbead_first,rsev,koffx,rcapmyp,l_break_sq,r_ring,wr,binding_rate_myo2,binding_rate_myp2,binding_rate_for,vpol,geometry);
% CAUTION: bmmat is not updated here. run bm_nbr again to update bmmat.
% % %     % record number of crosslinkers before breaking
% % %     nbefore = sum(xmat(:)) * .5;
    % indices of beads
    ibead = 1:size(rbead,2);
    % make it into vector form
    [xrow,xcol,xval] = find(xmat);
    % the default form of bpmat_trial, without any turnover, is equal to
    % bmmat
    bpmat_trial = bpmat;


        %% myosin turnover
        % number of myosin clusters to remove
        nmyoturn = koffmyo * size(rmyo,2) * dt;
        nmypturn = koffmyp * size(rmyp,2) * dt;
        % make it an integer
        nmyoturn = floor(nmyoturn) + (rand < nmyoturn - floor(nmyoturn));
        nmypturn = floor(nmypturn) + (rand < nmypturn - floor(nmypturn));
        % indices of myosin clusters to remove
        imyooff = randperm(size(rmyo,2), nmyoturn);
        % get rid of myo2 that fly too far
        temp = find(abs(rmyo(3,:)) > wr * 3);
        imyooff = unique([imyooff,temp]);
            % myo2 nodes take formins with them when they unbind
            temp = sum(fmmat(:,imyooff),2);
            ifiloff = find(temp)';
        imypoff = randperm(size(rmyp,2), nmypturn)
        % myp2 clusters that do not bind to any actin
        free_myp2 = sum(bpmat) ==0;
        free_myp2 = find(free_myp2);
        imypoff = union(free_myp2,imypoff)
        nmypfreeoff = numel(imypoff) - nmypturn;
        if ~isempty(rmyp)
            rmyp(:,imypoff) = [];
        end 
        % modify
        rmyo(:,imyooff) = [];
        bancm(imyooff) = [];
        fmmat(:,imyooff) = [];
        % sort ifiloff from largest to smallest
        ifiloff = sort(ifiloff,'descend');
        % modify
        for i = ifiloff
            rbead(:,ifor(i):ipt(i)) = [];
            ibead(ifor(i):ipt(i)) = [];
            bancf(ifor(i):ipt(i)) = [];
            bpmat_trial(ifor(i):ipt(i),:) = []; 
            fmmat(i,:) = [];
        end 
        nbead = length(ibead);
        % update xmat
            % delete entries in xmat that correspond to disappeared beads
            [~,xrow] = ismember(xrow,ibead);
            [~,xcol] = ismember(xcol,ibead);
            xval = xval(and(xrow~=0,xcol~=0));
            xrownew = xrow(and(xrow~=0,xcol~=0));
            xcol = xcol(and(xrow~=0,xcol~=0));
            xrow = xrownew;
        % update ifor and ipt
            ifor(~ismember(ifor,ibead)) = [];
            ipt(~ismember(ipt,ibead)) = [];
            [~,ifor] = ismember(ifor,ibead);
            [~,ipt] = ismember(ipt,ibead);
        dbead_first(ifiloff) = [];

        
        nmyo2_in = binding_rate_myo2 * 2 * pi * r_ring * dt;                % number of myosin clusters to come in
        nmyo2_in = floor(nmyo2_in) + (rand < nmyo2_in - floor(nmyo2_in));   % make it an integer
        for i = 1:nmyo2_in
            for ibind = 1:100
                theta = 2 * pi * rand;
                [x_in, y_in] = pol2cart(theta,r_ring-d_myo);
                z_in = wr * (rand - .5);
                % it shouldn't overlap with existing myo2
                temp = rmyo - repmat([x_in;y_in;z_in],1,size(rmyo,2));
                temp = sum(temp.*temp);
                if all(temp > .1*rexc_n*rexc_n)
                    break
                end
                if ibind == 100
                    %error('cannot insert myo2')
                    break
                end
            end
            rmyo = [rmyo,[x_in; y_in; z_in]];
            bancm = [bancm, true];
            fmmat = [fmmat, false(size(fmmat,1),1)];
        end
        
        %% formin turnover
        % indices of beads
        ibead = 1:size(rbead,2);
        % number of formin dimers to remove
        nforturn = kofffor * dt * length(ifor);
        % make it an integer
        nforturn = floor(nforturn) + (rand < (nforturn - floor(nforturn)));
        % indices of filament to be removed
        ifiloff = randperm(length(ifor),nforturn);
        % also remove filaments with abnormally long actin segments
        temp = circshift(rbead,[0,-1]) - rbead; % vector from this bead to the next bead
        seglensq = sum(temp.*temp);             % length of each segment
        seglensq(ipt) = 0;
        temp = find(seglensq > 4*dbead*dbead);
        disp(numel(temp))
        if ~isempty(temp)
            temp = temp(1);
            temp = find(ifor<=temp);
            temp = temp(end);
            ifiloff = [ifiloff, temp];
            ifiloff = unique(ifiloff);
        end
        % sort ifiloff from largest to smallest
        ifiloff = sort(ifiloff,'descend');        
        % modify
        for i = ifiloff
            rbead(:,ifor(i):ipt(i)) = [];
            ibead(ifor(i):ipt(i)) = [];
            bancf(ifor(i):ipt(i)) = [];
            bpmat_trial(ifor(i):ipt(i),:) = []; 
            fmmat(i,:) = [];
        end 
        nbead = length(ibead);
        % update xmat
            % delete entries in xmat that correspond to disappeared beads
            [~,xrow] = ismember(xrow,ibead);
            [~,xcol] = ismember(xcol,ibead);
            xval = xval(and(xrow~=0,xcol~=0));
            xrownew = xrow(and(xrow~=0,xcol~=0));
            xcol = xcol(and(xrow~=0,xcol~=0));
            xrow = xrownew;
        % update ifor and ipt
            ifor(~ismember(ifor,ibead)) = [];
            ipt(~ismember(ipt,ibead)) = [];
            [~,ifor] = ismember(ifor,ibead);
            [~,ipt] = ismember(ipt,ibead);
        dbead_first(ifiloff) = [];
        %% formin binding and nucleation
        nfor_in = (binding_rate_for/koffmyo)* 2 * pi * r_ring - length(ifor);   % number of formin dimers to nucleate
        nfor_in = floor(nfor_in) + (rand < (nfor_in - floor(nfor_in)));     % make it an integer
        
        if ~isempty(rmyo)
            for i = 1:nfor_in

%                 r_new_for = rmyo(:,myo_ind);    
                                                
                % now pick a random direction in 3D, see
                % http://mathworld.wolfram.com/SpherePointPicking.html, Eq. 16
                for inuc = 1:100
                    myo_ind = randi(size(rmyo,2));  % pick a myosin to bind to
                    temp = sum(fmmat,1);              % total number of formins on each myosin
                    x = rand;
                    %if temp(myo_ind)>3              % cannot bind to myosin already occupied by 4 formin dimers
                    if rand <= temp/2
                        continue
                    end
                    [thm, ~] = cart2pol(rmyo(1,myo_ind),rmyo(2,myo_ind));
                    zm = rmyo(3,myo_ind);
                    rho = r_ring - d_for;
                    [xf,yf,zf] = pol2cart(thm,rho,zm);
                    r_new_for = [xf; yf; zf];
                    x = randn;
                    y = randn;
                    z = randn;
                    vrand = [x;y;z];
                    vrand = vrand / norm(vrand);
                    if vrand(1:2)' * r_new_for(1:2) > 0
                        vrand = - vrand;     % cannot nucleate into the membrane       
                    end
                    % squared distance from this to all actin beads
                    temp = rbead;
                    temp(:,ifor) = 0;
                    temp = temp - repmat(r_new_for + 0.1 * dbead * vrand, 1, size(temp,2));
                    temp = sum(temp.*temp);
                    % if vrand is not within 2*rex_act of any bead, end
                    if all(temp > (1*rex_act)^2)
                        break
                    end
                    if inuc == 100
                        error('cannot insert actin')
                    end
                end
                r_new_act = r_new_for + 0.1 * dbead * vrand;    % position of the new actin bead

                rbead = [rbead, r_new_for, r_new_act];     
                ipt = [ipt, nbead + 2];
                ifor = [ifor, nbead + 1];
                nbead = nbead + 2;
                bancf = [bancf, true, false];
                bpmat_trial = [bpmat_trial; false(2,size(bpmat_trial,2))];
                dbead_first = [dbead_first, 0.1 * dbead];
                fmmat = [fmmat; false(1,size(fmmat,2))];
                fmmat(end,myo_ind) = true;
            end
        end
    %% cofilin severing
        % total length of actin filament
        ltot = dbead * (nbead - length(ifor));
        % number of severing events
        nsev = ltot * rsev * dt;
        % make it an integer
        nsev = floor(nsev) + (rand < nsev - floor(nsev));
        % cannot sever at formin beads or the bead after (leaving a bare formin)
        legalsite = 1:nbead;
        legalsite(or(ismember(legalsite,ifor),ismember(legalsite-1,ifor))) = [];
        % sites of severing
        nsev = min(nsev,length(legalsite));
        isev = legalsite(randperm(length(legalsite),nsev));
        % sort isev from largest to smallest
        isev = sort(isev,'descend');
 
        % modify
        ibead = 1:size(rbead,2);
        bfor = ismember(ibead,ifor); % boolean vector of formins
        bpt = ismember(ibead,ipt);
        for i = isev
            % the pointed end of the severed filament
            iptsev = find(bpt(i:end), 1) + i - 1; % note: ipt must be already sorted from smallest to largest
            rbead(:,i:iptsev) = [];
            bpmat_trial(i:iptsev,:) = [];   
            ibead(i:iptsev) = [];
            bancf(i:iptsev) = [];
            bfor(i:iptsev) = [];
            bpt(i:iptsev) = [];
            bpt(i-1) = true;
        end
        ifor = find(bfor);
        ipt = find(bpt);
        
        nbead = size(rbead,2);  
        % update xmat
        % delete entries in xmat that correspond to disappeared beads
        [~,xrow] = ismember(xrow,ibead);
        [~,xcol] = ismember(xcol,ibead);
        xval = xval(and(xrow~=0,xcol~=0));
        xrownew = xrow(and(xrow~=0,xcol~=0));
        xcol = xcol(and(xrow~=0,xcol~=0));
        xrow = xrownew;
%% myp2 binding            
        nmyp2_in = binding_rate_myp2 * 2 * pi * r_ring * dt;                           % number of myosin clusters to come in
        nmyp2_in = floor(nmyp2_in) + (rand < nmyp2_in - floor(nmyp2_in));   % make it an integer
        nmyp2_in = nmyp2_in + nmypfreeoff;
        if ~isempty(rbead)
            for i = 1:nmyp2_in
                % pick a random location in a cylinder
                if geometry == 'cyl'
                    for itrial = 1:1000
                        if itrial == 1000
                            %error('cannot insert myp2')
                            r_new = [];
                            break
                        end
                        xin = (rand*2-1) * (r_ring-rcapmyp);
                        yin = (rand*2-1) * (r_ring-rcapmyp);
                        zin = (rand*2-1) * (r_ring-rcapmyp);
                        rinsq = xin*xin + yin*yin;
                        if rinsq > (r_ring-rcapmyp)*(r_ring-rcapmyp)||abs(zin) > 3*wr
                            continue
                        end
                        temp = rbead;
                        temp(:,ifor) = 1e6;
                        temp(1,:) = temp(1,:) - xin;
                        temp(2,:) = temp(2,:) - yin;
                        temp(3,:) = temp(3,:) - zin;
                        temp = sum(temp.*temp);
                        if any(temp < rcapmyp.*rcapmyp)
                            r_new = [xin;yin;zin];
%                             disp(itrial)
                            break
                        end
                    end
                else
                    error('this geometry has not been coded')
                end
                rmyp = [rmyp,r_new];
                bpmat_trial = [bpmat_trial, false(size(bpmat_trial,1),1)];
            end
        end
    %% crosslinker breaking due to over extension
        cross_dx = rbead(:,xrow) - rbead(:,xcol);
        cross_dx2 = sum(cross_dx .* cross_dx);
        ind_del = cross_dx2 > l_break_sq;
        xrow(ind_del) = [];
        xcol(ind_del) = [];
        xval(ind_del) = [];
        xmat = sparse(xrow,xcol,xval,nbead,nbead);
    %% crosslinker turnover in itself
    if ~isempty(rbead)
        % number of xlinkers to go off
        nxturn = koffx * dt * .5 * full(sum(sum(xmat)));
        % make it an integer
        nxturn = floor(nxturn) + (rand < nxturn - floor(nxturn));
        % random indices to go off
        ixoff = randperm(sum(sum(xmat))/2,nxturn);
        % delete these xlinkers
        [xrow,xcol,xval] = find(triu(xmat));
        xval(ixoff) = false;
        xmat = sparse(xrow,xcol,xval,nbead,nbead);
        xmat = or(xmat,xmat');
        % number of crosslinkers after breaking
        nafter = sum(xmat(:)) * .5;
    end
        %% actin polymerization (add bead)
        [xrow,xcol,xval] = find(xmat);
        poly_ind = find(dbead_first > 1.1 * dbead);
        for ipoly = numel(poly_ind) : -1 : 1            % reverse order, convenient for insersion operations
            i = poly_ind(ipoly);
            dbead_first(i) = dbead_first(i) - dbead;
            r1 = rbead(:,ifor(i));
            r2 = rbead(:,ifor(i)+1);
            dr = r2 - r1;
            r_new = r2 - dbead * dr / norm(dr);     % position of the new bead to add
            rbead = [rbead(:,1:ifor(i)), r_new, rbead(:,ifor(i)+1:end)];
            bancf = [bancf(1:ifor(i)), false, bancf(ifor(i)+1:end)];
            bpmat_trial = [bpmat_trial(1:ifor(i),:); false(1,size(bpmat_trial,2)); bpmat_trial(ifor(i)+1:end,:)];
            
            nbead = nbead + 1;
            ifor = ifor + 1;
            ifor(1:i) = ifor(1:i) - 1;              % these two lines are equivalent to ifor(i+1:end) += 1, but avoid the endpoint problem
            ipt(i:end) = ipt(i:end) + 1;
            xrow = xrow + (xrow > ifor(i));         % to update xmat, it is sufficient to update xrow and xcol
            xcol = xcol + (xcol > ifor(i));
        end
        xmat = sparse(xrow,xcol,xval,nbead,nbead);
        dbead_first = dbead_first + vpol * dt;
        if isempty(rbead)
            rbead = double.empty(3,0);
        end
        if isempty(rmyo)
            rmyo = double.empty(3,0);
        end
        %% binding of crosslinkers
        % rules: density of xlinkers assumed to be constant in time: rhoain. total number of xlinkers: rhoain * r_ring. Add xlinkers such that the total number matches this value.
        nadd = round(rhoain*r_ring*2*pi) - nafter;
        % find all possible pairs within connection length of ain1
        possible_pairs = false(nbead);
        temp = 1:nbead;
        ind_not_for = setdiff(temp,ifor);
        temp = rangesearch(rbead',rbead',rxbind);
        for i = ind_not_for
            possible_pairs(i,temp{i}) = true;
            possible_pairs(i,1:i) = false;
        end
        possible_pairs(:,ifor) = false;
        possible_pairs(xmat)=false;
        % total number of possible binding sites
        nsite = sum(sum(possible_pairs));
        % there must be enough sites
        if nsite < nadd
            error('nsite < nadd')
        else
            selected = randsample(nsite,nadd);
        end
        % add to xmat
        [i,j] = find(possible_pairs);
        xrow = [xrow; i(selected); j(selected)];
        xcol = [xcol; j(selected); i(selected)];
        xval = [xval; true(nadd*2,1)];
        xmat = sparse(xrow,xcol,xval,nbead,nbead);
end
