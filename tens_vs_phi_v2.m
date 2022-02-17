% change from v1: include experimental values

clear
% preallocate
ilist = 1:10;
tlist = 1:33;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
        filename = strcat('dpfc1d_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        if exist(filename,'file')
                load(filename)
            tension_circ
%             tension_laplace_v1
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

tens_store_original = tens_store;
tens_store = tens_store_original * 1.9;

tens_mean = mean(tens_store);
tens_sd = std(tens_store);
%% plot tension versus phi
r_ring_t = 1.85 - 0.07*(tlist-10);
phi = acosd(r_ring_t / 1.85);
phi(tlist<10)=nan;
errorbar(real(phi),tens_mean,tens_sd,'-ok')
xlabel('\phi (deg)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
% axis([0,26,0,1600])
%% plot experimental values of tension versus phi
tens_exp = [425,535,555,690,605,800,770]; % wt
tens_exp_err = [133,270,390,334,183,268,198]; % wt
phi_exp = 10:10:70; % wt

tens_exp = [234, 453, 336, 292, 482, 526]; % del myp2
tens_exp_err = [270, 686, 277, 328, 226, 350]/2; % del myp2
phi_exp = 10:10:60; % del myp2

hold on 
errorbar(phi_exp,tens_exp,tens_exp_err,'bo');
hold off