% change from v1: include experimental values
% change from v2: use normalized ring length
% change from v3: use time

% clear
% preallocate
figure
ilist = 1:32;
tlist = 0:21;
% tlist = 0:60:2400;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for it = 1:numel(tlist)
        t = tlist(it);
%         filename = strcat('nosep3_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('wt2/fh2d9_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('slowp/slowp_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('wtfc4/wtfc4_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
%         filename = strcat('allmyo3/allmyo3_', num2str(ilist(itemp)), '_', num2str(t), 'sec.mat');
        filename = strcat('delmyp_fub/delmyp_fub_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        if exist(filename,'file')
                load(filename)
%             tension_circ
            tension_laplace_v1
            % averaging
            tens_store(itemp, it) = mean(tens);
        else
            disp('no file')
            break
        end
    end
end
% bad_cols = any(isnan(tens_store));
% tens_store(:,bad_cols) = [];
% tlist(bad_cols) = [];

% tens_store = tens_store * 1.5;

tens_mean = nanmean(tens_store, 1);
tens_sd = nanstd(tens_store, 0, 1);
%% plot tension versus phi
it = tlist;
errorbar(it(1:end-2),tens_mean(1:end-2),tens_sd(1:end-2),'ok', 'linewidth', 2)
xlabel('Time (min)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
axis([0,26,0,Inf])
%% plot experimental values of tension versus phi
tens_exp = [425,535,555,690,605,800,770];
tens_exp_err = [133,270,390,334,183,268,198];
phi_exp = 10:10:70;
t_exp = phi_exp/90 * 26;
hold on 
errorbar(t_exp,tens_exp,tens_exp_err,'bo', 'linewidth', 2);
set(gca,'Position',[0.2 0.2 0.7 0.7]);
set(gca, 'Linewidth', 4)
set(gca, 'FontSize', 30)
hold off