cols = {'#0072BD','#D95319','#000000','#7E2F8E','#77AC30'};
cnt = 1;
for k=0:10:40
%for k=0:5:20
	v = get_velocity_sliding_window('coopf',k,50,200,5,4:4:80);
	disp(mean(v))
	histogram(v*1000,'Normalization','probability','DisplayStyle','Stairs',...
			'EdgeColor',cols{cnt},'BinEdges',0:1:50,'LineWidth',2);
	cnt = cnt + 1;
	hold on;
end
[hleg,att] = legend('0','10','20','30','40');
%[hleg,att] = legend('0','5','10','15','20');
%title(hleg,'Window (s)')
title('Cdc12 node velocities (20 runs each)');
set(gca,'FontSize',14);
xlabel('Cdc12 Node velocity (nm/s)');
ylabel('prob. density');
