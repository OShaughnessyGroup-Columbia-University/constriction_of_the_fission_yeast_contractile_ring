%% decide what files are available
a = dir;
prefix = 'sizep5_';
iload = 64;
% axis_1 = [2,3,4,5]; 
% axis_2 = [1,3,5,7];
figure
timelist = (1:16);
for itime = 1:16
    time = timelist(itime);
    filename = strcat(prefix,num2str(iload),'_0min.mat');
    if ~exist(filename,'file')
        continue
    end
    load(filename,'rmyo','rmyp','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
    filename = strcat(prefix,num2str(iload),'_',num2str(time),'min.mat');
    if ~exist(filename,'file')
        continue
    end
    load(filename,'rmyo','rmyp','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
    subplot(4,4,itime)
%     fluo_myo_and_myp_script
    ringplot_toscale
    box_size = 2;
    axis([-box_size,box_size,-box_size,box_size'])
    title(['#',num2str(iload)])
%     plot_r_typical_v7
%     title([num2str(sep_myo_myp*1000,'%2.0f'), ' nm separation, ', num2str(act_spread*1000,'%2.0f'), ' nm spread'])
%     title([num2str(sep_myo_myp*1000,2), ' nm, ', num2str(bound_frac,2), ' bound'])
    view(0,90)
    axis off
end