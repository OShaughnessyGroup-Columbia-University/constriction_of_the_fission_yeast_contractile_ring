tens_100 = nan(7,100);
for iload = 1:100
    for jload = 2:7
        filename = ['drag2d2_', num2str(iload), '_', num2str(jload),'min.mat']
        load(filename)
        tension_circ
        tens_100(jload,iload) = mean(tens);
    end
end

tens_100 = tens_100(:,end:-1:1);
tens_100(tens_100 == 0) = nan;

plot(nanmean(tens_100))
xlabel('Number of steps in 2 s')
ylabel('Tension (pN)')
title('Average of 2 min to 7 min')

figure
x1 = tens_100(:,91:100);
x1 = nanmean(x1);
x2 = tens_100(:,6:15);
x2 = nanmean(x2);
boxplot([x2',x1'],'Labels',{'6 to 15','91 to 100'})
ylim([0,1000])
ylabel('Mean tension from 2 to 7 min (pN)')
[h,p] = ttest2(x1,x2);
title(['p = ', num2str(p), ', two-tailed t test'])

figure
x1 = tens_100(:,91:100);
x1 = x1(:);
x2 = tens_100(:,11);
boxplot(padcat(x2,x1),'Labels',{'#11','91 to 100'})
ylim([0,1000])
ylabel('Tension values every 1 min, from 2 to 7 min (pN)')
[h,p] = ttest2(x1,x2);
title(['p = ', num2str(p), ', two-tailed t test'])