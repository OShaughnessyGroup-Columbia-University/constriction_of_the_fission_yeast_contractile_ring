% change from v1: include experimental values
% change from v2: use normalized ring length
% change from v3: use time

% clear
% preallocate
% figure
% ilist = 1:24;
% tlist = 0:26;
ilist = 1:100;
tlist = 60:60:600;
% tlist = 1:26;
tens_store = nan(numel(ilist), numel(tlist));
tcirc_store = nan(numel(ilist), numel(tlist));
tf_store = nan(numel(ilist), numel(tlist));
ts_store = nan(numel(ilist), numel(tlist));
v_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    disp(['load run irun = ' num2str(ilist(itemp))])
    for it = 1:numel(tlist)
        disp(['t = ' num2str(tlist(it))])
%         clear vmyp
%         filename = strcat('nosep3_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('wt2/fh2d9_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('vconst2/vcont2_', num2str(ilist(itemp)), '_',...
%             num2str(t), 'sec.mat');
%         filename = strcat('delta_ain1p/delta-ain1p', num2str(ilist(itemp)), ...
%             '_', num2str(t), 'min.mat');
        filename = strcat('gm_2b/gm_2_', num2str(ilist(itemp)), '_', ...
            num2str(tlist(it)), 'sec.mat');
%         filename = strcat('gm_2/gm_2_', num2str(ilist(itemp)), '_', ...
%             num2str(tlist(it)), 'sec.mat');
%         filename = strcat('gm4/gm4_', num2str(ilist(itemp)), '_', ...
%             num2str(tlist(it)), 'sec.mat');
%         filename = strcat('gmdrop/gmdrop_', num2str(ilist(itemp)), '_', ...
%             num2str(tlist(it)), 'sec.mat');
%         filename = strcat('tent/tent_', num2str(ilist(itemp)), '_', ...
%             num2str(tlist(it)), 'sec.mat');
%         filename = strcat('noagm3/noagm3_', num2str(ilist(itemp)), '_', ...
%             num2str(tlist(it)), 'sec.mat');
        if exist(filename,'file')
            load(filename)
            tension_laplace_v1
            tens_store(itemp, it) = mean(tens);
%             tension_circ_v4
%             tcirc_store(itemp, it) = mean(tens);
            calc_sliding_fixed
            tf_store(itemp, it) = mean(tfixed);
            ts_store(itemp, it) = mean(tslide);
            v_store(itemp, it) = vave;
        else
            disp([filename ' does not exist'])
            break
        end
    end
end
% tens_store = tcirc_store;
% bad_cols = any(isnan(tens_store));
% tens_store(:,bad_cols) = [];
% tlist(bad_cols) = [];


tens_mean = nanmean(tens_store);
tens_sd = nanstd(tens_store);
%% plot tension versus phi
it = tlist;
figure
errorbar(it(1:end),tens_mean(1:end),tens_sd(1:end),'-o')
xlabel('Time (min)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)