temp = vbead(1:2,ifor);
vfor_model = sqrt(sum(temp.*temp)) * 1000; % unit: nm/s
temp = vmyo(1:2,:);
vmyo_model = sqrt(sum(temp.*temp)) * 1000; % unit: nm/s

vfor_mean_exp = 28;
vfor_sd_exp = 11;
vmyo_mean_exp = 22;
vmyo_sd_exp = 10;

errorbar_groups([mean(vfor_model),mean(vmyo_model);vfor_mean_exp,vmyo_mean_exp],[std(vfor_model),std(vmyo_model);vfor_sd_exp,vmyo_sd_exp], 'bar_names',{'Cdc12p','Myo2p'})