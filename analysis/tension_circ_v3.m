% change from original: use simplified version (segment tension)
% change from v2: consider the sign of segtens
%% Calculate tension
get_segtens                     % calculate the tension in each segment of actin filament (between beads)
rdiff = get_rdiff(rbead,1,ipt);
temp = circshift(rdiff,1,2);
segtens_s = -sum(temp.*segtens);
segtens_s(ifor+1) = segtens_s(ifor+1) .* dbead_first / 0.1;
tens = sum(segtens_s) * .1 / 2 / pi / r_ring