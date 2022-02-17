%% calculate steady state actin filament length probability vs filament length, given turnvoer parameters and lfvec,
% a vector containing filament lengths
function l_act_prob = l_act_dist(lfvec,rsev,vpol)
    kforoff = 0.023;
    temp = ((kforoff + lfvec * rsev)/vpol).*exp(-(lfvec.*(2*kforoff+lfvec*rsev))/2/vpol);
    l_act_prob = temp / sum(temp);
end
