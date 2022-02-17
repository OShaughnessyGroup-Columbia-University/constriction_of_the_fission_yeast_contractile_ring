clear
load tens_sizep4_lapv1
tens_store = nan(64,16);
for i = 1:64
    if numel(tens_cell{i}) >= 16
        tens_store(i,:) = tens_cell{i}(1:16);
    else
        tens_store(i,1:numel(tens_cell{i})) = tens_cell{i};
    end
end
means = 640; % experimental mean 640
sds = 290;   % experimental sd 290
exp_tens = [1118.403347,1194.508327,913.0030149,572.8077378,536.1879887,743.5950668,303.3453871,422.5820727,915.6224515,405.8880941,627.1723723,481.0612775,411.6186922,1019.993924,302.0408212,697.7888518,1100.081017,230.9960297,284.4688428,487.6357742,314.7616807,1022.710435,550.6567061,472.3519604,1093.54553,345.6881026,897.7164751,584.3531421,475.3916677,720.0088999,603.4689861];
for i = [1,2,3,4]
    ilist = (i*8-7):(i*8);
    temp = tens_store(ilist,:);
    means = [means, nanmean(temp(:))];
    sds = [sds, nanstd(temp(:))];
    [h,p]=ttest2(exp_tens,temp(:));
    p
end

bar(means)
hold on
errorbar(means,sds)
hold off
