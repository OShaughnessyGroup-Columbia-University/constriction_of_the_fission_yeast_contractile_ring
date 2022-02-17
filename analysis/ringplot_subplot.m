%% decide what files are available
a = dir;
prefix = 'fh2d2_220_';
suffix = 'min.mat';
i_six = 1:8
% figure
% for i = 1:8
%     iload = i_six(i);
%     filename = strcat(prefix,num2str(iload),'_',num2str(time),'min.mat');
%     if ~exist(filename,'file')
%         continue
%     end
%     load(filename,'rmyo','rmyp','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
%     subplot(2,4,i)
%     ringplot_pick
%     view(0,0)
% end
figure
for i = 1:8
    iload = i_six(i);
    filename = strcat(prefix,num2str(iload),suffix);
    if ~exist(filename,'file')
        continue
    end
    load(filename,'rmyo','rmyp','r_ring','fc','ipt','ifor','rbead','bancf','bancm','xmat')
    subplot(2,4,i)
%     fluo_myo_and_myp_script
    ringplot_toscale
%     ringplot_pick
    box_size = 2.9;
    axis([-box_size,box_size,-box_size,box_size'])
%     plot_r_typical
    title([num2str(iload-10),' min'])
    view(0,90)
end