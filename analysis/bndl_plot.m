%% preparing variables
tag = 'wtcon';
load([tag '_1_00sec.mat'])

bndl_store = [];
nx_store = [];
t_store = [];
w_store = [];
sep_store = [];

snap_time = 180;
nval3 = 1%size(X1, 3); % number of copies
nval2 = 1%size(X1, 2) % values of X1
nval1 = 1%size(X1, 1) % values of X2
nval = nval1*nval2
nmax = nval*12;
% ilist = 1:nval:nval2*nval;
ilist = 1:nval;
bndl_mat = nan(nval1, nval2);
width_mat = nan(nval1, nval2);
thick_mat = nan(nval1, nval2);
sep_mat = nan(nval1, nval2);
%% perform calculations
unit = 'sec';
for i0 = ilist
    disp(['******************** i0 = ' num2str(i0) ' ********************'])
    rHistogramWrap
    nx_store = [nx_store; nx_list];
    bndl_store = [bndl_store; bndl_act_frac];
    sep_store = [sep_store; myo_myp_sep*1e3];
    w_store = [w_store; width*1e3];
    t_store = [t_store; thick*1e3];
end

% nval = 1;
% nmax = 10;
% tlist = 30:30:1200;
% nt = length(tlist);
% bndl_mat = nan(nval, nt);
% width_mat = nan(nval, nt);
% thick_mat = nan(nval, nt);
% sep_mat = nan(nval, nt);
% i0 = 1;
% for snap_time = tlist
%     rHistogramWrap
%     bndl_store = [bndl_store; bndl_act_frac];
%     sep_store = [sep_store; myo_myp_sep*1e3];
%     w_store = [w_store; width*1e3];
%     t_store = [t_store; thick*1e3];
% end

width_mat(:) = nanmean(w_store, 2);
bndl_mat(:) = nanmean(bndl_store, 2);
bndl_sd = nanstd(bndl_store, 0, 2);
thick_mat(:) = nanmean(t_store, 2);
thick_sd = nanstd(t_store, 0, 2);
sep_mat(:) = nanmean(sep_store, 2);
sep_sd = nanstd(sep_store, 0, 2);

%% plotting the result
if(nval2 == 1 || nval1 == 1)
    figure;
    %bar3(thick_mat(2:end, :))
    bar(thick_mat)
    title('Bundle thickness')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
    
    % figure;
    % bar3(width_mat)
    % title('Bundle width')
    % set(gca, 'Linewidth', 3)
    % set(gca, 'FontSize', 24)
    % set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    % tit = get(gca,'Title');
    % view([-2 2 1])
    % xlim([0, nval2+1])
    % ylim([0, nval1+1])
    % saveas(gcf, ['./plots/' tit.String '.png'])
    % saveas(gcf, ['./plots/' tit.String '.fig'])
    
    figure;
    % bar3(bndl_mat(2:end, :))
    bar(100-bndl_mat)
    title('Whisker fraction')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
    
    figure;
    bar(sep_mat)
    title('Myo-Myp separation')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
else
    figure;
    %bar3(thick_mat(2:end, :))
    bar3(thick_mat)
    title('Bundle thickness')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    view([-2 2 1])
    xlim([0, nval2+1])
    ylim([0, nval1+1])
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
    
    % figure;
    % bar3(width_mat)
    % title('Bundle width')
    % set(gca, 'Linewidth', 3)
    % set(gca, 'FontSize', 24)
    % set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    % tit = get(gca,'Title');
    % view([-2 2 1])
    % xlim([0, nval2+1])
    % ylim([0, nval1+1])
    % saveas(gcf, ['./plots/' tit.String '.png'])
    % saveas(gcf, ['./plots/' tit.String '.fig'])
    
    figure;
    % bar3(bndl_mat(2:end, :))
    bar3(100-bndl_mat)
    title('Whisker fraction')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    view([-2 -2 1])
    xlim([0, nval2+1])
    ylim([0, nval1+1])
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
    
    figure;
    bar3(sep_mat)
    title('Myo-Myp separation')
    set(gca, 'Linewidth', 3)
    set(gca, 'FontSize', 24)
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    tit = get(gca,'Title');
    view([-2 2 1])
    xlim([0, nval2+1])
    ylim([0, nval1+1])
    saveas(gcf, ['./plots/' tit.String '.png'])
    saveas(gcf, ['./plots/' tit.String '.fig'])
end