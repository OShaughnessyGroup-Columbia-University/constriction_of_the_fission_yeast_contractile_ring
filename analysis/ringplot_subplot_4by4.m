%% decide what files are available
a = dir;
prefix = 'sizep4_';
time = 10;
axis_1 = [2,4,6,8]; 
axis_2 = [2,4,6,8];
% i_sixteen = find(ismember(X1,x1(axis_1))&ismember(X2,x2(axis_2)))
i_sixteen = [16,17,36,59];
figure
for i = 1:4
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
    subplot(1,4,i)
    fluo_myo_and_myp_script
%     ringplot_toscale
%     box_size = 1.9;
%     axis([-box_size,box_size,-box_size,box_size'])
%     plot_r_typical_v7
%     title(['#',num2str(iload)])
%     title([num2str(sep_myo_myp*1000,'%2.0f'), ' nm separation, ', num2str(act_spread*1000,'%2.0f'), ' nm spread'])
%     title([num2str(sep_myo_myp*1000,'%2.0f'), ' nm separation, ', num2str(bound_frac*100,'%2.0f'), '% bound'])
%     view(0,90)
%     axis off
end