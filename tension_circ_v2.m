% change from original: use simplified version (segment tension 
%% Calculate tension
get_segtens                     % calculate the tension in each segment of actin filament (between beads)
segtens_s = sqrt(sum(segtens.*segtens));
segtens_s(ifor+1) = segtens_s(ifor+1) .* dbead_first / 0.1;
tens = sum(segtens_s) * .1 / 2 / pi / r_ring