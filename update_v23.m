function [rbead,rmyo,rmyp,xmat,fmmat,ipt,ifor,bancf,bancm,dbead_first,bpmat_trial,t_rbead,t_rmyo,t_rmyp] =...
    update_v19(rbead,rmyo,rmyp,bpmat,xmat,rhoain,rxbind,rex_act,rexc_n,fmmat,bancf,bancm,ifor,ipt,kofffor,...
    koffmyo,koffmyp,dt,d_for,d_myo,dbead,dbead_first,rsev,koffx,rcapmyp,l_break_sq,r_ring,wr,...
    binding_rate_myo2,binding_rate_myp2,binding_rate_for,vpol,geometry,current_time,t_rbead,t_rmyo,t_rmyp,...
	alphafor,epsilonfor, qmax);
% CAUTION: bmmat is not updated here. run bm_nbr again to update bmmat.
% % %     % record number of crosslinkers before breaking
% % %     nbefore = sum(xmat(:)) * .5;

	if isempty(t_rmyo)
		t_rmyo = -10000*ones(1,size(rmyo,2));
    end
	
	if isempty(t_rbead)
		t_rbead = -10000*ones(1,size(rbead,2));
    end
	
	if isempty(t_rmyp)
		t_rmyp = -10000*ones(1,size(rmyp,2));
    end

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
        imypoff = randperm(size(rmyp,2), nmypturn);
        % myp2 clusters that do not bind to any actin
        free_myp2 = sum(bpmat) ==0;
        free_myp2 = find(free_myp2);
        imypoff = union(free_myp2,imypoff);
        nmypfreeoff = numel(imypoff) - nmypturn;
		
        % modify
	rmyo(:,imyooff) = [];
	t_rmyo(:,imyooff) = [];
		
	rmyp(:,imypoff) = [];
	t_rmyp(:,imypoff) = [];
		
        bancm(imyooff) = [];
        fmmat(:,imyooff) = [];
        % sort ifiloff from largest to smallest
        ifiloff = sort(ifiloff,'descend');
        % modify
        for i = ifiloff
            rbead(:,ifor(i):ipt(i)) = [];
	    t_rbead(ifor(i):ipt(i)) = [];
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
            qn = cart2pol(rmyo(1, :), rmyo(2, :));
            cts = histcounts(qn, -pi:0.1/r_ring:pi);
            empt = find(cts==0);
            if(~isempty(empt))
                theta = -pi + 0.1/r_ring*(empt(randperm(length(empt), 1))-0.5);
		[x_in, y_in] = pol2cart(theta,r_ring-d_myo);
		z_in = wr * (rand - .5);
            else
                ntrial = 100;
                for ibind = 1:ntrial
                    theta = 2 * pi * rand;
                    [x_in, y_in] = pol2cart(theta,r_ring-d_myo);
                    z_in = wr * (rand - .5);
                    % it shouldn't overlap with existing myo2
                    temp = rmyo - repmat([x_in;y_in;z_in],1,size(rmyo,2));
                    temp = sum(temp.*temp);
                    if all(temp > .1*rexc_n*rexc_n)
                        break
                    end
                    if ibind == ntrial
                        error('cannot insert myo2')
                    end
                end
            end
            rmyo = [rmyo,[x_in; y_in; z_in]];
	    t_rmyo = [t_rmyo,current_time + rand/10000];
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
        if ~isempty(temp)
            temp = temp(1);
            temp = find(ifor<=temp);
            temp = temp(end);
            ifiloff = [ifiloff, temp];
            ifiloff = unique(ifiloff);
        end
        % sort ifiloff from largest to smallest
        ifiloff = sort(ifiloff,'descend');        
        disp(['update removing ' num2str(numel(ifiloff)) ' overstretched fils'])
        % modify
        for i = ifiloff
            rbead(:,ifor(i):ipt(i)) = [];
			t_rbead(ifor(i):ipt(i)) = [];
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
		nfornode = sum(fmmat,1);  % total number of formins on each myosin
		
		% count number of each node and calculate probabilities of formin binding
		nxfor = histc(nfornode,[0 1 2]); 
		[p0,p1] = evolve_node_numbers(nxfor(1),nxfor(2),nxfor(3),...
						dt,koffmyo,alphafor,epsilonfor,false);
						% note: p0, p1 are meant to be applied sequentially.
						% i.e. look at 0 nodes and convert to 1 with prob. p0.
						% then, look at 1 nodes and convert to 2 with prob. p1.
		
		% construct array of indices of nodes which will receive formins
        addFormIndices = [];
		for count1=1:length(nfornode)
			if nfornode(count1) == 0 & rand() < p0
				addFormIndices = [addFormIndices count1];
				nfornode(count1) = nfornode(count1) + 1;
			end
			if nfornode(count1) == 1 & rand() < p1
				addFormIndices = [addFormIndices count1];
				nfornode(count1) = nfornode(count1) + 1;
			end
		end

        dbead0 = dbead;
        if ~isempty(rmyo)
            for myo_ind = addFormIndices
                for inuc = 1:400

                    %if(temp(myo_ind) > 0)
                    %    idx_for = find(fmmat(:, myo_ind), 1);
                    %    idx_bead = ifor(idx_for);
                    %    tfor = rbead(:, idx_bead+1) - rbead(:, idx_bead);
                    %    tring = [-rbead(2, idx_bead); rbead(1, idx_bead); 0];
                    %    polNew = sum(tfor.*tring); % correlated polarity
                    %else
                    %    polNew = sign(randn());
                    %end

                    [thm, ~] = cart2pol(rmyo(1,myo_ind),rmyo(2,myo_ind));
                    zm = rmyo(3,myo_ind);
                    rho = r_ring - d_for;
                    [xf,yf,zf] = pol2cart(thm,rho,zm);
                    r_new_for = [xf; yf; zf];
                    Nnewfil = 1;

                    r_new_act = zeros(3,Nnewfil);
                    % figure out angle that corresponds to 100 nm and polarity
                    actRodAng = (dbead/r_ring);
                    rhohat = [ xf; yf; 0]/rho;
                    phihat = [-yf; xf; 0]/rho;
                    zhat   = [  0;  0; 1];
		    cqmax = cos(qmax);
                    for k=1:Nnewfil
		        polNew = sign(randn()); % random polarity
                        cpsi = rand()*(1 - cqmax) + cqmax;
                        psi = acos(cpsi);
                        q = rand()*pi;
                        tNew = polNew*cpsi*phihat + ...
                               sin(psi)*(cos(q)*zhat - sin(q)*rhohat);
                        r_new_act(:,k) = r_new_for + tNew*dbead;
                    end

                    % squared distance from this to all actin beads
                    temp = rbead;
                    temp(:,ifor) = 0;
                    temp = temp - repmat(r_new_act, 1, size(temp,2));
                    temp = sum(temp.*temp);
                    % if vrand is not within 2*rex_act of any bead, end
                    if all(temp > (1*rex_act)^2)
                        break
                    end
                    if inuc == 400
                        'cannot insert actin'
			continue;
                    end
                end
		        rbead = [rbead, r_new_for, r_new_act];
		        t_rbead = [t_rbead, current_time + rand(1,Nnewfil + 1)/10000];				
		        ipt = [ipt, nbead + Nnewfil + 1];
		        ifor = [ifor, nbead + 1];
		        nbead = nbead + Nnewfil + 1;
		        bancf = [bancf, true, false(1,Nnewfil)];
		        bpmat_trial = [bpmat_trial; false(Nnewfil + 1,size(bpmat_trial,2))];
		        dbead_first = [dbead_first, dbead0];
		        fmmat = [fmmat; false(1,size(fmmat,2))];
		        fmmat(end,myo_ind) = true;
            end
        end
    %% cofilin severing
        % old severing scheme
        % % total length of all actin filaments
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

        % % biased severing towards pointed end
        % isev = zeros(1, nsev);
        % for i = 1:nsev
        %     for itrial = 1:100
        %         ibd = randi(nbead);
        %         idxfor = find(ifor <= ibd, 1, 'last');
        %         if rand() < (ibd - ifor(idxfor) - 1)*0.01
        %             isev(i) = ibd;
        %             break;
        %         else
        %             continue;
        %         end
        %     end 
        % end

        
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
			t_rbead(i:iptsev) = [];
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
                    for itrial = 1:3000
                        if itrial == 3000
                            error('cannot insert myp2')
                        end
                        xin = (rand*2-1) * (r_ring-rcapmyp);
                        yin = (rand*2-1) * (r_ring-rcapmyp);
                        zin = (rand*2-1) * 3*wr;
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
				t_rmyp = [t_rmyp,current_time + rand/10000];
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
			t_rbead = [t_rbead(1:ifor(i)), current_time + rand/10000, t_rbead(ifor(i)+1:end)];
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
			t_rbead =  double.empty(1,0);
        end
        if isempty(rmyo)
            rmyo = double.empty(3,0);
			t_rmyo = double.empty(1,0);
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
            disp('warning: nsite < nadd')
        else
            selected = randsample(nsite,nadd);
	    % add to xmat
	    [i,j] = find(possible_pairs);
	    xrow = [xrow; i(selected); j(selected)];
	    xcol = [xcol; j(selected); i(selected)];
	    xval = [xval; true(nadd*2,1)];
	    xmat = sparse(xrow,xcol,xval,nbead,nbead);
        end
end
