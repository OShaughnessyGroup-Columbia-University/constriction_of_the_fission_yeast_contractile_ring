% change from v1: include experimental values
% change from v2: use normalized ring length
% change from v3: use time

% clear
% preallocate
figure
ilist = 1:64;
dt = 1;
tlist = 2;
% unit = 'sec';
unit = 'min';
% tlist = 1:26;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
%         filename = strcat('nosep3_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('wt2/fh2d9_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        filename = strcat('wtfc3/wtfc3_', num2str(ilist(itemp)), '_', num2str(t), unit, '.mat')
%         filename = strcat('matur/matur_', num2str(ilist(itemp)), '_', num2str(t), unit, '.mat')
%         filename = strcat('allmyo/allmyo_', num2str(ilist(itemp)), '_', num2str(t), 'sec.mat')
        if exist(filename,'file')
                load(filename)
            tension_laplace_v1
%             tension_circ_v4
            % averaging
%             tens_store(itemp, t) = mean(tens);
            tens_store(itemp, t/dt) = mean(tens);
        else
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
t = tlist;
errorbar(t(1:end-2),tens_mean(1:end-2),tens_sd(1:end-2),'-o')
xlabel('Time (min)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
% axis([0,26,0,Inf])
%% plot experimental values of tension versus phi
% tens_exp = [425,535,555,690,605,800,770];
% tens_exp_err = [133,270,390,334,183,268,198];
% phi_exp = 10:10:70;
% t_exp = phi_exp/90 * 26;
% hold on 
% errorbar(t_exp,tens_exp,tens_exp_err,'bo');
% hold off