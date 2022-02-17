% average tension over different runs and over time.

tArray = 50:5:250;
iArray = 1:4:80;
%iArray = 1:1:20;
tensArray  = nan(length(iArray),length(tArray));

for t1=1:length(tArray)
for k=1:length(iArray)
	f1 = ['coopf_' num2str(iArray(k)) '_' num2str(tArray(t1)) 'sec.mat'];
	if isfile(f1)
		load(f1,'fanc','rbead','rmyo','rmyp','r_ring');
		tension_laplace_v1;
		tensArray(k,t1) = tens;
	end
end
end

plot(tArray,nanmean(tensArray),'bo')