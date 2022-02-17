function l_t = l_act_t(lfvec,rsev,vpol,t)
    l_t = sum(lfvec .* (1-exp(-(rsev*lfvec*t))) ./ (rsev * lfvec *t) .* l_act_dist(lfvec,rsev,vpol));
end