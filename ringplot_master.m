%% decide what files are available
clear
a = dir;
prefix = 'fcold_3_';
suffix = 'min.mat';

% load file
for i = 1:2:100
    if exist(strcat(prefix,num2str(i),suffix),'file')
        load(strcat(prefix,num2str(i),suffix),'rmyo','rmyp','r_ring','fc','nbead','ipt','ifor','rbead','bancf','bancm','xmat')
    figure
    ringplot_toscale
    view(0,90)
    end
end