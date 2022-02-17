%% change from original: use different ga for Myo2p and Myp2p

function amat = get_a_v2(rbead, rmyo, rmyp, bmmat, bpmat, gppf, rdiff, dt, kcap, kcap_p, ifor, ipt, bmmatpt, bpmatpt, ga_rat)
% get the (sparse) A matrix (see notes "det ring new scheme")
%% preallocation
    % Estimate: there are nmyo myosin clusters, each interacts with ~20
    % filaments, so ~ 20*nmyo pairs. Each pair gives 18 non-zero entries in
    % the A matrix, so 360 * nmyo non-zero elements. There are nbead beads,
    % each gives 18 non-zero entries due to artificial drag, so 18 * nbead
    % non-zero elements. Use a larger number for safety.
nb = size(rbead,2);
nm = size(rmyo,2) + size(rmyp,2);
nmyo2 = size(rmyo,2);
nmyp2 = size(rmyp,2);
%% f-v relation 
% bead-list and myosin-list of interacting bead-myosin pairs
[blist, mlist] = find([bmmat,bpmat]);
lenblist = length(blist);
%% respective ga for myo2 and myp2
ga_myo2 = ga_rat * dt * kcap;
ga_myp2 = ga_rat * dt * kcap_p;
	%% BLK1
	ib = blist;
	im = mlist + nb;
	g = gppf(mlist);
    g(mlist<=nmyo2) = g(mlist<=nmyo2) - ga_myo2;
    g(mlist>nmyo2) = g(mlist>nmyo2) - ga_myp2;
	rd = rdiff(:,ib);

	nn = zeros(9,length(ib));

	for i = 1:length(ib)
		t1 = rd(:,i) * rd(:,i)';
		nn(:,i) = t1(:);
    end
    
    if isempty(g)
        gmat = [];
    else
        gmat = nn .* g(ones(1,9),:);
    end
	aval1 = gmat(:)';
	aval2 = gmat(:)';
	
	arow1 = zeros(1,9*lenblist);
	acol1 = zeros(1,9*lenblist);
	
	arow2 = zeros(1,9*lenblist);
	acol2 = zeros(1,9*lenblist);
	
	arow1(1:9:9*lenblist) = 1 + 3 * (ib-1);
	arow1(2:9:9*lenblist) = 2 + 3 * (ib-1);
	arow1(3:9:9*lenblist) = 3 + 3 * (ib-1);
	arow1(4:9:9*lenblist) = 1 + 3 * (ib-1);
	arow1(5:9:9*lenblist) = 2 + 3 * (ib-1);
	arow1(6:9:9*lenblist) = 3 + 3 * (ib-1);
	arow1(7:9:9*lenblist) = 1 + 3 * (ib-1);
	arow1(8:9:9*lenblist) = 2 + 3 * (ib-1);
	arow1(9:9:9*lenblist) = 3 + 3 * (ib-1);
	
	arow2(1:9:9*lenblist) = 1 + 3 * (im-1);
	arow2(2:9:9*lenblist) = 2 + 3 * (im-1);
	arow2(3:9:9*lenblist) = 3 + 3 * (im-1);
	arow2(4:9:9*lenblist) = 1 + 3 * (im-1);
	arow2(5:9:9*lenblist) = 2 + 3 * (im-1);
	arow2(6:9:9*lenblist) = 3 + 3 * (im-1);
	arow2(7:9:9*lenblist) = 1 + 3 * (im-1);
	arow2(8:9:9*lenblist) = 2 + 3 * (im-1);
	arow2(9:9:9*lenblist) = 3 + 3 * (im-1);
	
	acol1(1:9:9*lenblist) = 1 + 3 * (im-1);
	acol1(2:9:9*lenblist) = 1 + 3 * (im-1);
	acol1(3:9:9*lenblist) = 1 + 3 * (im-1);
	acol1(4:9:9*lenblist) = 2 + 3 * (im-1);
	acol1(5:9:9*lenblist) = 2 + 3 * (im-1);
	acol1(6:9:9*lenblist) = 2 + 3 * (im-1);
	acol1(7:9:9*lenblist) = 3 + 3 * (im-1);
	acol1(8:9:9*lenblist) = 3 + 3 * (im-1);
	acol1(9:9:9*lenblist) = 3 + 3 * (im-1);
	
	acol2(1:9:9*lenblist) = 1 + 3 * (ib-1);
	acol2(2:9:9*lenblist) = 1 + 3 * (ib-1);
	acol2(3:9:9*lenblist) = 1 + 3 * (ib-1);
	acol2(4:9:9*lenblist) = 2 + 3 * (ib-1);
	acol2(5:9:9*lenblist) = 2 + 3 * (ib-1);
	acol2(6:9:9*lenblist) = 2 + 3 * (ib-1);
	acol2(7:9:9*lenblist) = 3 + 3 * (ib-1);
	acol2(8:9:9*lenblist) = 3 + 3 * (ib-1);
	acol2(9:9:9*lenblist) = 3 + 3 * (ib-1);

	%% BLK 1 is the vectorized version of following jos code
	%{
	for i = 1:length(blist)
		ib = blist(i);
		im = mlist(i) + size(rbead,2);
		g = gppf(mlist(i)) - ga; % gppf means gamma_pull per filament. It is a vector with nmyo elements.
		rd = rdiff(:,ib); % it is a column-vector
		% get the local gamma (3 X 3 matrix). It will be put into the sparse matrix A. 
		% for the equation, see notes "effective drag force"
		gmat = g * (rd * rd'); % note that gmat is symmetric
	%     [grow, gcol, gval] = find(gmat);
		% put gmat into the sparse matrix A 
		arow(indvec) = [1 2 3 1 2 3 1 2 3] + 3 * (ib-1);
		acol(indvec) = [1 1 1 2 2 2 3 3 3] + 3 * (im-1);
		aval(indvec) = gmat(:);
		indvec = indvec + 9;
		arow(indvec) = [1 2 3 1 2 3 1 2 3] + 3 * (im-1);
		acol(indvec) = [1 1 1 2 2 2 3 3 3] + 3 * (ib-1);
		aval(indvec) = gmat(:);
		indvec = indvec + 9;
		% update ind
		ind = ind + 18;
	end
	%}
	
	%% artificial drag between myosin and actin
	%% BLK 2
	
	[bptlist, mptlist] = find([bmmatpt,bpmatpt]);
	lenbptlist = length(bptlist);
	
	ib = bptlist;
	im = mptlist + nb;
	
	garow = zeros(1,6*lenbptlist);
	gacol = zeros(1,6*lenbptlist);
	gaval = zeros(1,6*lenbptlist);
	
	garow(1:6:6*lenbptlist) = 1 + 3 * (ib-1);
	garow(2:6:6*lenbptlist) = 2 + 3 * (ib-1);
	garow(3:6:6*lenbptlist) = 3 + 3 * (ib-1);
	garow(4:6:6*lenbptlist) = 1 + 3 * (im-1);
	garow(5:6:6*lenbptlist) = 2 + 3 * (im-1);
	garow(6:6:6*lenbptlist) = 3 + 3 * (im-1);
	
	gacol(1:6:6*lenbptlist) = 1 + 3 * (im-1);
	gacol(2:6:6*lenbptlist) = 2 + 3 * (im-1);
	gacol(3:6:6*lenbptlist) = 3 + 3 * (im-1);
	gacol(4:6:6*lenbptlist) = 1 + 3 * (ib-1);
	gacol(5:6:6*lenbptlist) = 2 + 3 * (ib-1);
	gacol(6:6:6*lenbptlist) = 3 + 3 * (ib-1);
	
	gaval(1:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
	gaval(2:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
    gaval(3:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
    gaval(4:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
    gaval(5:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
    gaval(6:6:6*lenbptlist) = ga_myo2 * (mptlist <= nmyo2) + ga_myp2 * (mptlist > nmyo2);
	%% BLK2 is the vectorized version of jos code below	
	%{
	garow = zeros(1,720*size(rmyo,2));
	gacol = garow;
	gaval = garow;
	ind2 = 1;
	gmat = ga * eye(3);
	[grow, gcol, gval] = find(gmat);
	[bptlist, mptlist] = find(bmmatpt);
	for i = 1:length(bptlist)
		ib = bptlist(i);
		im = mptlist(i) + size(rbead,2);
		n = length(grow);
		% put gmat into the sparse matrix ga
		garow(ind2:ind2+n-1) = grow' + 3 * (ib-1);
		gacol(ind2:ind2+n-1) = gcol' + 3 * (im-1);
		gaval(ind2:ind2+n-1) = gval';
		garow(ind2+n:ind2+2*n-1) = grow' + 3 * (im-1);
		gacol(ind2+n:ind2+2*n-1) = gcol' + 3 * (ib-1);
		gaval(ind2+n:ind2+2*n-1) = gval';
		% update ind
		ind2 = ind2 + 2*n;
	end
	%}
	
	%% artifical drag between adjacent actin beads
	%% if this bead is not a formin, or the next bead after a formin,
	%% or the pointed end
	
	%% BLK 3
	arow3 = zeros(1,3*nb);
	acol3 = zeros(1,3*nb);
	aval3 = zeros(1,3*nb);
	
	arow3(1:3:3*nb) = 1 + 3*((1:nb)-1);
	arow3(2:3:3*nb) = 2 + 3*((1:nb)-1);
	arow3(3:3:3*nb) = 3 + 3*((1:nb)-1);
	
	acol3(1:3:3*nb) = 1 + 3*((1:nb)-2);
	acol3(2:3:3*nb) = 2 + 3*((1:nb)-2);
	acol3(3:3:3*nb) = 3 + 3*((1:nb)-2);
	
	aval3 = repmat(ga_myo2,1,3*nb);
	
	aval3(3*(ifor-1)+1) = 0;
	aval3(3*(ifor-1)+2) = 0;
	aval3(3*(ifor-1)+3) = 0;
	
	aval3(3*(ifor+1-1)+1) = 0;
	aval3(3*(ifor+1-1)+2) = 0;
	aval3(3*(ifor+1-1)+3) = 0;
	
	
	arow4 = zeros(1,3*nb);
	acol4 = zeros(1,3*nb);
	aval4 = zeros(1,3*nb);
	
	arow4(1:3:3*nb) = 1 + 3*((1:nb)-1);
	arow4(2:3:3*nb) = 2 + 3*((1:nb)-1);
	arow4(3:3:3*nb) = 3 + 3*((1:nb)-1);
	
	acol4(1:3:3*nb) = 1 + 3*((1:nb));
	acol4(2:3:3*nb) = 2 + 3*((1:nb));
	acol4(3:3:3*nb) = 3 + 3*((1:nb));
	
	aval4 = repmat(ga_myo2,1,3*nb);
	
	aval4(3*(ifor-1)+1) = 0;
	aval4(3*(ifor-1)+2) = 0;
	aval4(3*(ifor-1)+3) = 0;
	
	aval4(3*(ipt-1)+1) = 0;
	aval4(3*(ipt-1)+2) = 0;
	aval4(3*(ipt-1)+3) = 0;	
	
	%% END OF BLK3
	
	%% BLK 3 is a vectorized version of the following jos code
	%{
	for i = 1:nb
		% if this bead is not a formin or the next bead after a formin
		if ~ or(any(ifor == i), any(ifor+1 == i))
			% artificial drag between this bead and the previous bead
			arow(ind:ind+2) = (1:3) + 3*(i-1);
			acol(ind:ind+2) = (1:3) + 3*(i-2);
			aval(ind:ind+2) = repmat(ga,1,3);
	%         arow(ind+3:ind+5) = (1:3) + 3*(i-2);
	%         acol(ind+3:ind+5) = (1:3) + 3*(i-1);
	%         aval(ind+3:ind+5) = repmat(ga,1,3);
			ind = ind + 6;
		end
		% if this bead is not a pointed end or a formin
		if ~or(any(ipt == i), any(ifor == i))
			% artifical drag between this bead and the next bead
			arow(ind:ind+2) = (1:3) + 3*(i-1);
			acol(ind:ind+2) = (1:3) + 3*(i);
			aval(ind:ind+2) = repmat(ga,1,3);
	%         arow(ind+3:ind+5) = (1:3) + 3*(i);
	%         acol(ind+3:ind+5) = (1:3) + 3*(i-1);
	%         aval(ind+3:ind+5) = repmat(ga,1,3);
			ind = ind + 6;
		end
	end
	%}

	%aval(arow==0) = [];
	%arow(arow==0) = [];
	%acol(acol==0) = [];
	%amat = sparse(arow, acol, aval, 3*(size(rbead,2)+size(rmyo,2)), 3*(size(rbead,2)+size(rmyo,2)), ind-1);
	%gaval(garow==0) = [];
	%garow(garow==0) = [];
	%gacol(gacol==0) = [];
	%gamat = sparse(garow, gacol, gaval, 3*(size(rbead,2)+size(rmyo,2)), 3*(size(rbead,2)+size(rmyo,2)), ind2-1);
	%amat = amat + gamat;
	%{
	DEBUG
	'row1'
	find(arow1 == 0)
	'row2'
	find(arow2 == 0)
	'row3'
	find(arow3 == 0)
	'row4'
	find(arow4 == 0)
	'garow'
	find(garow == 0)
	
	'col1'
	find(acol1 == 0)
	'col2'
	find(acol2 == 0)
	'col3'
	find(acol3 == 0)
	'col4'
	find(acol4 == 0)
	'gacol'
	find(gacol == 0)
	%}
	%% the 4:end for arow/col/val3 is because of faulty indexing
	% in acol for the first bead.
	
	combrow = [arow1 arow2 arow3(4:end) arow4 garow];
	combcol = [acol1 acol2 acol3(4:end) acol4 gacol];
	combval = [aval1 aval2 aval3(4:end) aval4 gaval];
	% DEBUG
	%max(combrow)
	%max(combcol)
	%3*(nb+nm)
	amat = sparse(combrow, combcol, combval, 3*(nb+nm), 3*(nb+nm));
end
