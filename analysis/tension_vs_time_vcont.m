% change from v1: include experimental values
% change from v2: use normalized ring length
% change from v3: use time

% clear
% preallocate
figure
hold on
for i0 = 1:6
    ilist = i0:6:60;
    tlist = 0:10:1500;
    % tlist = 1:26;
    tens_store = nan(numel(ilist), numel(tlist));
    for itemp = 1:numel(ilist)
        for it = 1:numel(tlist)
            %         filename = strcat('nosep3_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
            %         filename = strcat('wt2/fh2d9_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
            filename = strcat('vcont2/vcont2_', num2str(ilist(itemp)), '_',...
                num2str(tlist(it)), 'sec.mat');
            %         filename = strcat('delta_ain1p/delta-ain1p', num2str(ilist(itemp)), ...
            %             '_', num2str(t), 'min.mat');
            %         filename = strcat('ptanch/ptanch_', num2str(ilist(itemp)), '_', ...
            %             num2str(tlist(it)), 'sec.mat');
            if exist(filename,'file')
                load(filename)
                % tension_laplace_v1
                tension_circ
                % averaging
                %             tens_store(itemp, 1+t) = mean(tens);
                %             tens_store_lap(itemp, it) = mean(tens);
                tens_store(itemp, it) = mean(tens);
            else
                disp([filename ' does not exist'])
                break
            end
        end
    end
    % bad_cols = any(isnan(tens_store));
    % tens_store(:,bad_cols) = [];
    % tlist(bad_cols) = [];
    
    % tens_store = tens_store * 1.5;
    
    tens_mean = nanmean(tens_store);
    tens_sd = nanstd(tens_store);
    %% plot tension versus phi
    it = tlist;
    errorbar(it(1:end),tens_mean(1:end),tens_sd(1:end),'-o')
end
xlabel('Time (min)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)