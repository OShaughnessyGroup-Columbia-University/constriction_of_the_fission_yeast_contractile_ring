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
for i = [1,2,3,4]
    ilist = (i*8-7):(i*8);
    temp = tens_store(ilist,:);
    means = [means, nanmean(temp(:))];
    sds = [sds, nanstd(temp(:))];
end

bar(means)
hold on
errorbar(means,sds)
hold off
ttest2