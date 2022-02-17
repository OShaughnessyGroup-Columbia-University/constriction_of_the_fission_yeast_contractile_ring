%% decide what files are available
a = dir;
prefix = 'delmyp6_';
time = 9;
i_sixteen = 1:64;
figure
for i = 1:64
    iload = i_sixteen(i);
    filename = strcat(prefix,num2str(iload),'_00min.mat');
    if ~exist(filename,'file')
        continue
    end
    load(filename,'rmyo','rmyp','rcapmyp','rcapmyo_short','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
    filename = strcat(prefix,num2str(iload),'_',num2str(time),'min.mat');
    if ~exist(filename,'file')
        continue
    end
    load(filename,'rmyo','rmyp','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
    subplot(8,8,i)
%     fluo_myo_and_myp_script
%     ringplot_toscale
%     box_size = 0.8;
%     axis([-box_size,box_size,-box_size,box_size'])
    plot_r_typical
    title([num2str(sep_myo_myp*1000,'%2.0f'), ' nm separation, ', num2str(act_spread*1000,'%2.0f'), ' nm spread'])
%     view(0,90)
    axis off
end