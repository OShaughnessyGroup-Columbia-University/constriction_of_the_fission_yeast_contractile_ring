function [] = run_and_plot(filename)
    loc_min = strfind(filename,'min');
    suffix = filename(loc_min:end);
    temp = strfind(filename,'_');
    loc_last_underscore = temp(end);
    prefix = filename(1:loc_last_underscore);
    para_file = [prefix,'0',suffix];
    load(para_file)
    load(filename)
    
    nstep = 10;
    dt = 0.02;
    rk1_v13
    ringplot_toscale
end