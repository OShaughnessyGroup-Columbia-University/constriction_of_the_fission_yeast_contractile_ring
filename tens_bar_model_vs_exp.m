temp = tens_store(:,2:end);
tens_mean_model = mean(temp(:))
tens_sd_model = std(temp(:));
temp = [320,484,478,290,558,310,51,80,774,851,322,153,650,401,135,565,354,338];
tens_mean_exp = mean(temp);
tens_sd_exp = std(temp);
bar([tens_mean_model,tens_mean_exp],.5,'c')
hold on
errorbar([tens_mean_model,tens_mean_exp],[tens_sd_model,tens_sd_exp],'k.')
hold off
set(gca,'XTickLabel',{'Model','Experiment'},'fontsize',20)
ylabel('Tension (pN)','FontSize',20)