% change from original: no force between formin and the first actin bead
% change from v1: no excluded volume force on the first actin bead whatsoever
% change from v2: force always perpendicular to local tangent vector
function F = act_excl_force_v1(rbead,rdiff,kex,rex,rcut,ifor)

	rbeadTr = rbead';
	nbead = size(rbead,2);

	% method 2
	% find all beads within distance rcut of each other
	[nList,dd] = rangesearch(rbeadTr,rbeadTr,rcut);
	
	F = zeros(3, nbead);
	
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
		
		% calculate mutual normal vector of neighbors except the self-term
		rr = cross(repmat(rdiff(:, k),[1, N-1]), rdiff(:, nListVec(2:N)));
    for nl = find(sum(rr.^2) < 1e-6)
        rr(:, nl) = rbead(:, nListVec(1+nl)) - rbead(:, k);
    end
		rhat = rr ./ repmat(sqrt(sum(rr.*rr)), [3, 1]);

    delr = repmat(rbead(:, k),[1, N-1]) - rbead(:, nListVec(2:N));
    delr = delr ./ repmat(sqrt(sum(delr.*delr)), [3, 1]);
		
		% calculate the force, excluding the self term
		F1 = kex * exp(-(ddVec(2:N).^2)./rex^2) .* sum(rhat.*delr);
		F(:, k) = rhat * F1';
	end
	
	F(:,ifor) = 0;
  F(:,ifor+1) = 0;
end
