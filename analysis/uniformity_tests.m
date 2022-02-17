

%% preparing variables
tag = 'gm6';
it = 180;
load([tag '_1_00sec.mat'])

dphi = 2*pi/30;
p = dphi/2/pi;

nval3 = 10% size(X1, 3); % number of copies
nval2 = 1 %size(X1, 2) % values of X1
nval1 = size(X1, 1) % values of X2
nval = nval1*nval2
nmax = nval*10

boxo_store = nan(nval1, nval2, nval3);
boxp_store = nan(nval1, nval2, nval3);

for irun = 1:nmax
    try
        load([ pfx '_' num2str(irun) '_' num2str(it) 'sec.mat']);
    catch
        continue
    end
    phimyo = cart2pol(rmyo(1, :), rmyo(2, :));
    phimyp = cart2pol(rmyp(1, :), rmyp(2, :));
    %% boxing method
    rho = histcounts(phimyo, -pi:dphi:pi);
    % disp('Boxing method myo2 (>1 means clustered)')
    boxo_store(irun) = std(rho)/sqrt((1-p)*length(phimyo)*p);
    rho = histcounts(phimyp, -pi:dphi:pi);
    % disp('Boxing method myp2 (>1 means clustered)')
    boxp_store(irun) = std(rho)/sqrt((1-p)*length(phimyp)*p);
    % figure
    % hold on
    
    % histogram(phimyp, -pi:dphi:pi)
    % histogram(phimyo, -pi:dphi:pi)
end
boxo_mat = nanmean(boxo_store, 3);
boxp_mat = nanmean(boxp_store, 3);

%% entropy method
% rho = rho/length(rmyo);
% disp('Entropy relative to uniform dist (<0 means more clustered)')
% S = sum(-rho.*log(rho)) - log(1/p)
% 
% %% nearest neighbor method
% [idx, dd ] = knnsearch(phimyo', phimyo', 'k', 2);
% disp('Mean nearest neighbor dist/expected (< 1 means clustered)')
% dnn = mean(dd(:, 2))/(pi/length(rmyo))