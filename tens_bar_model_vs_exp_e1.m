temp = tens_store(2:end);
tens_mean_model = nanmean(temp(:))
tens_sd_model = nanstd(temp(:));
temp = [124,127,62,264,155,171,153,188,147];
tens_mean_exp = mean(temp);
tens_sd_exp = std(temp);
bar([tens_mean_model,tens_mean_exp],.5,'c')
hold on
errorbar([tens_mean_model,tens_mean_exp],[tens_sd_model,tens_sd_exp],'k.')
hold off
set(gca,'XTickLabel',{'Model','Experiment'},'fontsize',20)
ylabel('Tension (pN)','FontSize',20)