% change from v1: include experimental values
% change from v2: use normalized ring length

clear
% preallocate
ilist = 1:10;
tlist = 1:35;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
        filename = strcat('fcnew_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
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

tens_store = tens_store * 1.5;

tens_mean = mean(tens_store);
tens_sd = std(tens_store);
%% plot tension versus phi
r_ring_t = 1.85 - 0.07*(tlist-10);
x = r_ring_t / 1.85;
errorbar(x,tens_mean,tens_sd,'-ok')
set(gca, 'xdir','reverse')
xlabel('\phi (deg)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
axis([0,1,0,Inf])
%% plot experimental values of tension versus phi
tens_exp = [425,535,555,690,605,800,770];
tens_exp_err = [133,270,390,334,183,268,198];
phi_exp = 10:10:70;
x_exp = cosd(phi_exp);
hold on 
errorbar(x_exp,tens_exp,tens_exp_err,'bo');
hold off