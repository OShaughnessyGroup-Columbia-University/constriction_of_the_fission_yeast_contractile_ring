% change from v1: bug fixed -- crosslinkers did not bind the ring in v1.

function [rbead,rmyo,rmyp,xmat,fmmat,ipt,ifor,bancf,bancm,dbead_first,bpmat_trial] = update_v2(rbead,rmyo,rmyp,bpmat,xmat,rhoain,rxbind,fmmat,bancf,bancm,ifor,ipt,kofffor,...
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
        nmypturn = koffmyp * size(rmyp,2)/2 * dt;
        % make it an integer
        nmyoturn = floor(nmyoturn) + (rand < nmyoturn - floor(nmyoturn));
        nmypturn = floor(nmypturn) + (rand < nmypturn - floor(nmypturn));
        % indices of myosin clusters to remove
        imyooff = randperm(size(rmyo,2), nmyoturn);
            % myo2 nodes take formins with them when they unbind
            temp = sum(fmmat(:,imyooff),2);
            ifiloff = find(temp)';
        imypoff = randperm(size(rmyp,2), nmypturn);
        % myp2 clusters that do not bind to any actin
        free_myp2 = sum(bpmat) ==0;
        free_myp2 = find(free_myp2);
        imypoff = union(free_myp2,imypoff);
        nmypfreeoff = numel(imypoff) - nmypturn;
        % modify
        rmyo(:,imyooff) = [];
        rmyp(:,imypoff) = [];
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
            theta = 2 * pi * rand;
            if geometry == 'cyl'
                [x_in, y_in] = pol2cart(theta,r_ring-d_myo);
                z_in = wr * (rand - .5);
                rmyo = [rmyo,[x_in; y_in; z_in]];
            end
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
        nfor_in = binding_rate_for * 2 * pi * r_ring * dt;                     % number of formin dimers to nucleate
        nfor_in = floor(nfor_in) + (rand < (nfor_in - floor(nfor_in)));     % make it an integer
        
        if ~isempty(rmyo)
            for i = 1:nfor_in
                myo_ind = randi(size(rmyo,2));  % pick a myosin to bind to
                temp = sum(fmmat,1);              % total number of formins on each myosin
                if temp(myo_ind)>3              % cannot bind to myosin already occupied by 4 formin dimers
                    continue
                end
                [thm, ~] = cart2pol(rmyo(1,myo_ind),rmyo(2,myo_ind));
                zm = rmyo(3,myo_ind);
                rho = r_ring - d_for;
                [xf,yf,zf] = pol2cart(thm,rho,zm);
                r_new_for = [xf; yf; zf];
%                 r_new_for = rmyo(:,myo_ind);    
                                                
                % now pick a random direction in 3D, see
                % http://mathworld.wolfram.com/SpherePointPicking.html, Eq. 16
                x = randn;
                y = randn;
                z = randn;
                vrand = [x;y;z];
                vrand = vrand / norm(vrand);
                if vrand(1:2)' * vrand(1:2) > r_ring*r_ring
                    vrand = - vrand;     % cannot nucleate into the membrane       
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
            
        nmyp2_in = binding_rate_myp2 * 2 * pi * r_ring * dt;                           % number of myosin clusters to come in
        nmyp2_in = floor(nmyp2_in) + (rand < nmyp2_in - floor(nmyp2_in));   % make it an integer
        nmyp2_in = nmyp2_in + nmypfreeoff;
        if ~isempty(rbead)
            for i = 1:nmyp2_in
                ind = 1:size(rbead,2);
                ind(ifor) = [];                     % indices for actin beads (non-formin)
                ind = ind(randi(numel(ind)));  % pick one at random
                for itrial = 1:1e3
                    x = randn;      % pick a random direction in 3D
                    y = randn;
                    z = randn;
                    vrand = [x;y;z];
                    vrand = vrand / norm(vrand);
                    % position of new myp2
                    r_new = rbead(:,ind) + rcapmyp * vrand;
                    if r_new(1:2)' * r_new(1:2) < r_ring*r_ring   % stop trying if some r_new that's within the membrane is found
                        break
                    end
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
        nadd = round(rhoain*r_ring) - nafter;
        % find all possible pairs within connection length of ain1
        possible_pairs = false(nbead);
        if nadd > 0
            for i = 1:nbead-1
                if any(ifor==i)
                    continue
                end
                for j = i+1:nbead
                    if any(ifor==j)
                        continue
                    end
                    if norm(rbead(:,i) - rbead(:,j)) < rxbind && ~xmat(i,j)
                        possible_pairs(i,j) = true;
                    end
                end
            end 
        end
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