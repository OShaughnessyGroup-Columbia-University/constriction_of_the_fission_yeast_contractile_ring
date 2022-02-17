% change from v1: include experimental values
% change from v2: include tension_laplace_v1 as another option. include del
% myp2 as another option
% change from v5: generate a 4 by 4 plot
% change from v6: use tens_cell in fh2d2_all_512_tens_vs_phi.mat as input
% v7 to v10 skipped
% change from v11: add best fit line of experimental data to the plot
% change from v12: calculate MSE from simulated tension to best fit line of
% experiment
% change from v13: bin model predictions, use chi-squared instead of MSE

%% control the weird color of errorbars
opengl hardwarebasic
%% decide what files are available
a = dir;
% time = 2;
axis_1 = 1:2:7; % row indices
axis_2 = 1:2:7; % column indices
[A1, A2] = meshgrid(axis_1,axis_2);
tlist = 0:16;

load('tens_myo2e1_circ.mat')
load('harvey_wt_tens_vs_phi.mat')

for i = 1:1
    disp(i)
    subplot(4,4,i)
%     ilist = find(ismember(X1,x1(A1(i)))&ismember(X2,x2(A2(i))));    
    ilist = 1:64;    
    tens_store = nan(8,length(tlist));
    for j = 1:length(ilist)
        if length(tlist) > length(tens_cell{ilist(j)})
            tens_store(j,1:length(tens_cell{ilist(j)})) = tens_cell{ilist(j)};
        else
            tens_store(j,:) = tens_cell{ilist(j)}(1:length(tlist));
        end
    end
    
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
    hold on
    % plot(real(phi),tens_store,'-ok');
%     xlabel('\phi (deg)','FontSize',18)
%     ylabel('Tension (pN)','FontSize',18)
    % axis([0,26,0,1600])
    %% plot experimental values of tension versus phi
    tens_exp = [425,535,555,690,605,800,770]; % wt
    tens_exp_err = [133,270,390,334,183,268,198]; % wt
    phi_exp = 10:10:70; % wt
    
%     tens_exp = [234, 453, 336, 292, 482, 526]; % del myp2
%     tens_exp_err = [270, 686, 277, 328, 226, 350]/2; % del myp2
%     phi_exp = 10:10:60; % del myp2

    % plot linear fit
    temp = [0,80];
    plot(temp,b(1)+b(2)*temp,'b');
    
    % errorbar plot of harvey's data
    errorbar(phi_exp,tens_exp,tens_exp_err,'bo');
    hold off
%     axis([0,90,0,1500])    
    
    %% bin model predictions
    phi_binned = discretize(phi(~isnan(phi)),[0,15,25,35,45,55,65,90]);
    sigma = tens_exp_err(phi_binned);
    x = tens_store(:,~isnan(phi));
    mu = tens_exp(phi_binned);
%     temp = (x - repmat(mu,8,1)).^2 ./ repmat(sigma.^2,8,1);
    temp = (x - repmat(mu,8,1)).^2 ./ nanstd(x).^2;
    chi_sq = nansum(temp(:));
    text(60,1800,['\chi^{2} = ',num2str(chi_sq,'%3.0f')])
    axis([-Inf,Inf,0,2000])
end