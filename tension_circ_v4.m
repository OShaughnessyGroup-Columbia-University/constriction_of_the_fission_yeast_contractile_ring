% change from original: use simplified version (segment tension)
% change from v2: consider the sign of segtens
% change from v3: consider the angular positions of beads
%% Calculate tension
get_segtens                     % calculate the tension in each segment of actin filament (between beads)
rdiff = get_rdiff(rbead,1,ipt);
temp = circshift(rdiff,1,2);
segtens_s = -sum(temp.*segtens);
segtens_s(ifor+1) = segtens_s(ifor+1) .* dbead_first / 0.1;
[th,~] = cart2pol(rbead(1,:),rbead(2,:));
dth = abs(- circshift(th,1,2) + th);
dth(dth > pi) = 2*pi-dth(dth > pi);
tens = sum(segtens_s.*dth) / 2 / pi;