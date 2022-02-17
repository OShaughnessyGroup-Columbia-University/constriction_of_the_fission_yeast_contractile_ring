% change from v1: include experimental values
% change from v2: include tension_laplace_v1 as another option. include del
% myp2 as another option
% change from v5: generate a 4 by 4 plot

%% decide what files are available
a = dir;
prefix = 'delmyp6_';
% time = 2;
% axis_1 = [1,2,3,4]; % row indices
% axis_2 = [1,3,6,8]; % column indices
% [A1, A2] = meshgrid(axis_1,axis_2);
tlist = 0:16;

for i = 1
    disp(i)
    subplot(4,4,i)
    ilist = 1:64;
%     ilist = find(ismember(X1,x1(A1(i)))&ismember(X2,x2(A2(i))));
    tens_store = nan(numel(ilist), numel(tlist));
    for itemp = 1:numel(ilist)
        for t = tlist
            filename = strcat(prefix, num2str(ilist(itemp)), '_00min.mat');
            if exist(filename,'file')
                load(filename)
            end
            filename = strcat(prefix, num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
            if exist(filename,'file')
                    load(filename)
%                 tension_circ
                tension_laplace_v1
                % averaging
                tens_store(itemp,t+1) = mean(tens);
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

    if size(tens_store,1) == 1
        tens_mean = tens_store;
        tens_sd = zeros(size(tens_store));
    else
        tens_mean = nanmean(tens_store);
        tens_sd = nanstd(tens_store);
    end

    %% plot tension versus phi
    r_ring_t = 1.85 - 0.07*(tlist-0);
    phi = acosd(r_ring_t / 1.85);
    phi(tlist<0)=nan;
    errorbar(real(phi),tens_mean,tens_sd,'-ok')
    % plot(real(phi),tens_store,'-ok');
    xlabel('\phi (deg)','FontSize',18)
    ylabel('Tension (pN)','FontSize',18)
    % axis([0,26,0,1600])
    %% plot experimental values of tension versus phi
%     tens_exp = [425,535,555,690,605,800,770]; % wt
%     tens_exp_err = [133,270,390,334,183,268,198]; % wt
%     phi_exp = 10:10:70; % wt

    tens_exp = [234, 453, 336, 292, 482, 526]; % del myp2
    tens_exp_err = [270, 686, 277, 328, 226, 350]/2; % del myp2
    phi_exp = 10:10:60; % del myp2

    hold on 
    errorbar(phi_exp,tens_exp,tens_exp_err,'bo');
    hold off
%     axis([0,90,0,1500])
end
temp = tens_store(:,2:17);
nanmean(temp(:))
nanstd(temp(:))