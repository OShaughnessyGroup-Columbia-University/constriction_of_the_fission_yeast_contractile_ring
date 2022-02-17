% change from original: no force between formin and the first actin bead
% change from v1: no excluded volume force on the first actin bead whatsoever
function F = act_excl_force_v1(rbead,kex,rex,rcut,ifor)

	rbeadTr = rbead';
	nbead = size(rbead,2);

	% method 1 -- ignore
	%{
	[nList,dd] = knnsearch(rbeadTr,rbeadTr,'K',N);
		
	F1 = kex * exp(-(dd.^2)./rex^2);
	F2 = zeros(nbead,3);
	
	for k = 1:nbead
		rr = repmat(rbeadTr(k,:),[N-1 1]) - rbeadTr(nList(2:N),:);
		rhat = rr./repmat(sqrt(sum(rr.*rr,2)),[1 3]);
		F2(k,:) = sum(F1(k,2:N) * rhat,1);
	end
	%}
	
	% method 2
	
	% find all beads within distance rcut of each other
	[nList,dd] = rangesearch(rbeadTr,rbeadTr,rcut);
	
	F2 = zeros(nbead,3);
	
	% loop through each bead
	for k = 1:nbead
        % skip if this bead is a formin
        if any(ifor==k) || any(ifor+1)==k
            continue
        end
        
        % locate neighborlist and distance list
		nListVec = nList{k};
		ddVec = dd{k};
        
		% find number of neighbors
		N = length(nListVec);
		
		% calculate rhat to all neighbors except the self-term
		rr = repmat(rbeadTr(k,:),[N-1 1]) - rbeadTr(nListVec(2:N),:);
		rhat = rr./repmat(sqrt(sum(rr.*rr,2)),[1 3]);
		
		% calculate the force, excluding the self term
		F1 = kex * exp(-(ddVec(2:N).^2)./rex^2);
		F2(k,:) = sum(F1 * rhat,1);
	end
	
	F = F2';
	F(:,ifor) = 0;
    F(:,ifor+1) = 0;
end
