% change from v1: include experimental values
% change from v2: include tension_laplace_v1 as another option. include del
% myp2 as another option
% change from v5: generate a 4 by 4 plot
% change from v6: no plot, save tension as a cell structure

%% decide what files are available
a = dir;
time = 2;
axis_1 = 1:8; % row indices
axis_2 = 1:8; % column indices
[A1, A2] = meshgrid(axis_1,axis_2);
tlist = 1:33;

tens_mean_store = nan(512,33);
tens_sd_store = nan(512,33);
for i = 1:64
    disp(i)
    
    ilist = find(ismember(X1,x1(A1(i)))&ismember(X2,x2(A2(i))));
    tlist = 1:33;
    tens_store = nan(numel(ilist), numel(tlist));
    for itemp = 1:numel(ilist)
        for t = tlist
            filename = strcat('fh2d2_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
            if exist(filename,'file')
                    load(filename)
                tension_circ
    %             tension_laplace_v1
                % averaging
                tens_store(itemp,t) = mean(tens);
            else
    %             break
            end
        end
    end
    tens_mean(i,:) = mean(tens_store);
    tens_sd(i,:) = std(tens_store);
end