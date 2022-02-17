% change from v1: include experimental values
% change from v2: include tension_laplace_v1 as another option. include del
% myp2 as another option
% change from v5: generate a 4 by 4 plot
% change from v6: no plot, save each run to a separate .mat file (total 512
% files)
% v7 and v8 skipped

%% decide what files are available
a = dir;
% time = 2;
axis_1 = [1,3,6,8]; % row indices
axis_2 = [1,3,6,8]; % column indices
[A1, A2] = meshgrid(axis_1,axis_2);
tlist = 1:27;

for ilist = 1:512
    ilist
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
    % bad_cols = any(isnan(tens_store));
    % tens_store(:,bad_cols) = [];
    % tlist(bad_cols) = [];

    % tens_store_original = tens_store;
    % tens_store = tens_store_original;

    save(['tens_',num2str(ilist),'.mat'],'tens_store')
end