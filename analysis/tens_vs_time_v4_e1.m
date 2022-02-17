% change from v1: include experimental values
% change from v2: use normalized ring length
% change from v3: use time

clear
figure
% preallocate
ilist = 12;
tlist = 1:35;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
        filename = strcat('e1_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
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

tens_store = tens_store;

tens_mean = mean(tens_store);
tens_sd = std(tens_store);
%% plot tension versus phi
t = tlist-10;
length(t)
length(tens_mean)
errorbar(t(1:end-2),tens_mean(1:end-2),tens_sd(1:end-2),'-ok')
xlabel('Time (min)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
axis([0,26,0,Inf])
%% plot experimental values of tension versus phi
tens_exp = [425,535,555,690,605,800,770];
tens_exp_err = [133,270,390,334,183,268,198];
phi_exp = 10:10:70;
t_exp = phi_exp/90 * 26;
hold on 
errorbar(t_exp,tens_exp,tens_exp_err,'bo');
hold off