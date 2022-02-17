clear
% preallocate
ilist = 1:10;
tlist = 1:26;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
        filename = strcat('fcold_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        if exist(filename,'file')
                load(filename)
            tension_circ
            % averaging
            tens_store(itemp,t) = mean(tens);
        else
            break
        end
    end
end
bad_cols = any(isnan(tens_store));
tens_store(:,bad_cols) = [];
tlist(bad_cols) = [];
tens_mean = mean(tens_store);
tens_sd = std(tens_store);
% plot tension versus time
errorbar(tlist,tens_mean,tens_sd)