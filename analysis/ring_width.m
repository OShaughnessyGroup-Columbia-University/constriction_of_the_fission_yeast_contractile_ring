workingDir = '.';
mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','sizep5', 10)]);
% mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','myo2e2', 10)]);
% mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','delmyp7', 10)]);
nmax = length(mat_files);
% nmax = 1;
dmyo_ave = zeros(1,nmax);
dmyp_ave = zeros(1,nmax);
dact_ave = zeros(1,nmax);
myo_w = zeros(1,nmax);
myo_t = zeros(1,nmax);
myp_w = zeros(1,nmax);
myp_t = zeros(1,nmax);
act_w = zeros(1,nmax);
act_t = zeros(1,nmax);
rwhisk = zeros(1,nmax);
std75 = 1.15035; % convert from stdev to 75%ile
% dcut = d_myo+2*rcapmyp;
dcut = 0.3;
% 
% figure
% hold all
%% 
for i=1:nmax
    load(mat_files(i).name);
    
    rbead(:,ifor) = []; %#ok<*SAGROW>
    ract = sqrt(rbead(1,:).^2 + rbead(2,:).^2);
    dact = r_ring - ract;
    dmean = median(dact);
    whiskers = (dact > dmean + 2*rcapmyp);
    rwhisk(i) = length(find(whiskers))/length(dact);
    dact(whiskers) = [];
    rbead(:,whiskers) = [];
    zmean = median(rbead(3,:));
    whiskers = ( abs(rbead(3,:) -zmean) > 2*rcapmyp);
    rwhisk(i) = rwhisk(i) + length(find(whiskers))/length(dact);
    dact(whiskers) = [];
    rbead(:,whiskers) = [];
    
    rmyo_r = sqrt(rmyo(1,:).^2 + rmyo(2,:).^2);
    rmyp_r = sqrt(rmyp(1,:).^2 + rmyp(2,:).^2);
    whiskers = (r_ring - rmyp_r > dmean + 2*rcapmyp);
    rmyp_r(whiskers) = [];
    rmyp(:, whiskers) = [];
    whiskers = (abs(rmyo(3,:) - zmean) > 2*rcapmyp);
    rmyp_r(whiskers) = [];
    rmyp(:, whiskers) = [];
    
    dmyo_ave(i) = median(r_ring - rmyo_r);
    dmyp_ave(i) = median(r_ring - rmyp_r);
    dact_ave(i) = median(dact);
    
    zmean = median(rbead(3,:));
    zact = rbead(3,:);
    nonWhiskerZ = (zact-zmean).^2 <= (2*rcapmyp)^2;
    zact(~nonWhiskerZ) = [];
    % Using quantiles
    dact_thick(i) = quantile(dact,0.95) - quantile(dact, 0.05);
    dact_width(i) = quantile(zact,0.95) - quantile(zact, 0.05);
    % Using variance -- for comparison with mcdonald
    zmyo = rmyo(3,:);
    nonWhiskerZ = (zmyo-zmean).^2 <= (2*rcapmyp)^2;
    zmyo(~nonWhiskerZ) = [];
    myo_w(i) = sqrt(var(zmyo) + rcapmyo_short^2);
    zmyp = rmyp(3,:);
    nonWhiskerZ = (zmyp-zmean).^2 <= (2*rcapmyp)^2;
    zmyp(~nonWhiskerZ) = [];
    myp_w(i) = sqrt(var(zmyp) + rcapmyp^2);
    myo_t(i) = sqrt(var(rmyo_r) + rcapmyo_short^2);
    myp_t(i) = sqrt(var(rmyp_r) + rcapmyp^2);
end
%%
act_thick = mean(dact_thick)
act_width = mean(dact_width)
rw_ave = mean(rwhisk)

myo_tave = mean(sqrt(myo_t.^2))*std75*2
myo_wave = mean(sqrt(myo_w.^2))*std75*2
myp_tave = mean(sqrt(myp_t.^2))*std75*2
myp_wave = mean(sqrt(myp_w.^2))*std75*2

sep = mean(dmyp_ave - dmyo_ave)