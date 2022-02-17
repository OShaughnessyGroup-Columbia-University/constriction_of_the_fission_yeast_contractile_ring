%% change from original: use different ga for Myo2p and Myp2p that are calcualted within the program

function m = get_m_v1 (dt,kcap,kcap_p,ifor,ipt,rbead,rmyo,rmyp,bancf,bancm,bmmat,bpmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt,bpmatpt)

ga_myo2 = 2 * 0.1 * kcap;
ga_myp2 = 2 * 0.1 * kcap_p;
% ga_myo2 = 2 * dt * kcap;
% ga_myp2 = 2 * dt * kcap_p;

% will be run repeatedly
% calculates the w matrix (see sec 2.5 in witkin et al)
% first w of each individual particle is calculated, then concatenated in
% block-diagonal form.
    nbead = size(rbead,2);
    % preallocate
        % total number of particles times 3
		nb=size(rbead,2);
		nm=size(rmyo,2);
        np=size(rmyp,2);
        if nb + nm + np == 0
            m = [];
            return
        end
    % for beads
    if ~isempty(rbead)
        totdrag_vec = ga_myo2 * sum(bmmatpt,2)+ ga_myp2 * sum(bpmatpt,2);
    else
        totdrag_vec = [];
    end
    
	
	% calculate nn for each bead
    nn = zeros(9,nb);
	
	for i = 1:nb
		t1 = rdiff(:,i) * rdiff(:,i)';        
        nn(:,i) = t1(:);        
    end
	
	
	%% artificial drag between adjacent actin beads on each filament
    if ~isempty(rbead)
        totdrag = 2 * ones(1,nb);
        totdrag(ifor) = totdrag(ifor) - 2;
        totdrag(ifor+1) = totdrag(ifor+1) - 1;
        totdrag(ipt) = totdrag(ipt) - 1;
    else
        totdrag = [];
    end
	
	% above lines implement
	% totdrag(i) = 2 - 2 * any(ifor == i) - any(ifor+1 == i) - any(ipt == i);
	
	%% zeroth order drag
    if ~isempty(rbead)
        kb = gb * ones(1,nb); % cytosolic drag
        kb(bancf) = kb(bancf) + gf; % membrane drag
    else
        kb = [];
    end
	kb = kb + ga_myo2 * totdrag + totdrag_vec'; 
		% first term: artificial drag between adjacent actin beads on each filament
		%% second term: artificial drag between myosin and actin beads

	% for myo2
	km = gmsol * ones(1,nm) +  ga_myo2 * sum(bmmatpt,1); % solution drag
    %% second term is artificial drag between myosin and actin beads
	km(bancm) = km(bancm) + gm; %if the ith myosin is anchored, membrane drag
	kp = gmsol * ones(1,np) + ga_myp2 * sum(bpmatpt,1);
    
	t1 = eye(3);
    if isempty(rbead)
        gbeadfv = [];
    end
	t2 = (gbeadfv-totdrag_vec)';
	A = bsxfun(@(x,y) y*x,t1(:),[kb km kp]);
	
	%size(t2)
	%size(nn)
	%size(t2(ones(1,9),:))
	if ~isempty(rbead)
        A(:,1:nb) = A(:,1:nb) + nn .* t2(ones(1,9),:);
    end
	
	%{
	A

A =

     1     2     3
     4     5     6
     7     8     9

>> B

B =

     1     4    10

   res =  A .* B(ones(1, 3),:);
   
     1     8    30
     4    20    60
     7    32    90
	%}
	
	%{
	for i = 1:nb
		t1 = kb(i) * eye(3);
		t2 = (gbeadfv(i)-totdrag_vec(i)) * nn(:,:,i);
        A{i} = sparse(t1 + t2);
			%% second term is effective drag from f-v relation
    end
	%}
	
    for i = 1:size(rmyo,2)
        % sum of nn for all beads interacting with the ith myosin
        nni = sum(nn(:,bmmat(:,i)),2); 
        % effective drag from f-v relation
        m = (gppf(i)-ga_myo2) * nni;
        % store
        A(:,nb+i) = A(:,nb+i) + m;
    end
		
    for i = 1:size(rmyp,2)
        % sum of nn for all beads interacting with the ith myosin
        nni = sum(nn(:,bpmat(:,i)),2); 
        % effective drag from f-v relation
        m = (gppf(nm+i)-ga_myp2) * nni;
        % store
        A(:,nb+nm+i) = A(:,nb+nm+i) + m;
    end
    
    % basically, you want cols to go 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 ...
	% and rows to go 1 2 3 1 2 3 1 2 3 4 5 6 4 5 6 4 5 6 ...
	% the below statements are just fancy ways of doing that
	cols = bsxfun(@(x,y) y*x,ones(1,3),(1:3*(nb+nm+np))')';
	
	t1 = bsxfun(@(x,y) y*x,ones(1,3),(1:3:3*(nb+nm+np))')';
	t2 = t1(:)';
    rows = [t2;t2+1;t2+2];
	
	m = sparse(rows(:)',cols(:)',A(:),3*(nb+nm+np),3*(nb+nm+np));
end