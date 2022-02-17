rHistogramWrap
x1 = 6:4:18;
x2 = 2:4:14;
[fun_myo, fun_myp] = meshgrid(x1,x2);

rbndl_2d= zeros(4,4);
rbndl_2d_std= zeros(4,4);
width_2d= zeros(4,4);
width_2d_std = zeros(4,4);
thick_2d= zeros(4,4);
thick_2d_std= zeros(4,4);
bndl_frac_2d= zeros(4,4);
bndl_frac_2d_std= zeros(4,4);
sep_2d= zeros(4,4);
sep_2d_std= zeros(4,4);
indx  = zeros(4,4);
for i = 1:16
    indx(i) = i;
    rbndl_2d(i) = mean(rbndl(i:16:end));
    thick_2d(i) = mean(thick(i:16:end));
    width_2d(i) = mean(width(i:16:end));
    bndl_frac_2d(i) = mean(bndl_act_frac(i:16:end));
    sep_2d(i) = mean(myo_myp_sep(i:16:end));
    rbndl_2d_std(i) = std(rbndl(i:16:end));
    thick_2d_std(i) = std(thick(i:16:end));
    width_2d_std(i) = std(width(i:16:end));
    bndl_frac_2d_std(i) = std(bndl_act_frac(i:16:end));
    sep_2d_std(i) = std(myo_myp_sep(i:16:end));
end

figure
for i = 1:16
    load(['bind_2da/bind_2da_' num2str(i) '_600sec.mat']);
    subplot(4,4,i)
    plot_cross_section_26oct18_fudge
end
