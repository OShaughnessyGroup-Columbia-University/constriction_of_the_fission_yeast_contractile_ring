clear
% preallocate
ilist = 1:64;
tlist = 1:26;
% tlist = 30:30:600;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
%         filename = strcat('g_anch_scan/g_anch_scan_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('v_sept_scan_fast_print/v_sept_scan_',
%         num2str(ilist(itemp)), '_', num2str(t), 'sec.mat');
        filename = strcat('wtfc3/wtfc3_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        if exist(filename,'file')
            load(filename)
            tension_laplace_v1
%             tension_circ_v4
            % averaging
            tens_store(itemp,t) = mean(tens);
        else
            break
        end
    end
end
tens_store
bad_cols = any(isnan(tens_store))
tens_store(:,bad_cols) = [];
tlist(bad_cols) = [];
tens_mean = mean(tens_store);
tens_sd = std(tens_store);
% plot tension versus time
figure
mean(tens_mean)
errorbar(tlist,tens_mean,tens_sd)