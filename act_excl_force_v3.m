% change from original: no force between formin and the first actin bead
% change from v1: no excluded volume force on the first actin bead whatsoever
% change from v2: switch to rod-based excluded volume method
function F = act_excl_force_v3(rbead,kex,rex,rcut,ifor)

	  rbeadTr = rbead';
	  nbead = size(rbead,2);

	  % find all beads within distance rcut of each other
	  [nList,dd] = rangesearch(rbeadTr,rbeadTr,rcut);
	  F = zeros(3, nbead);
	  
	  % loop through each bead
    for k = 1:nbead
        % skip if this segment contains a formin
        if any(ifor==k) || any(ifor+1==k)
          continue
        end
        
        % locate neighbor list and distance list
%         good = (nList{k} > k+1 );
        good = (nList{k} > k+1 & (ismember(nList{k}, nList{k-1}) ...
            | ismember(nList{k}-1, nList{k})) & ~ismember(nList{k}, ifor) ...
            & ~ismember(nList{k}, ifor+1));
        nListVec = nList{k}(good);
        ddVec = dd{k}(good);
        
        %loop over other segments
        for l = 1:length(nListVec)
          % s is normalized distance from brbded end bead to pt of closes approach
          % F1 is a 1x3 vector
          [F1, s, t] = getSegForce(k, nListVec(l), rbead, kex, rex);
          
          F(:, k) = F(:, k) + F1*s;
          F(:, k-1) = F(:, k-1) + F1*(1-s);
          F(:, nListVec(l)) = F(:, nListVec(l)) - F1*t;
          F(:, nListVec(l)-1) = F(:, nListVec(l)-1) - F1*(1-t);
        end
    end
	  %F = F2';
	  F(:,ifor) = 0;
    %F(:,ifor+1) = 0;
end
