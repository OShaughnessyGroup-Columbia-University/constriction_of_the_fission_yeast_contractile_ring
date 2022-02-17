%% decide what files are available
a = dir;
prefix = 'fcold_';
suffix = 'min.mat';

for jj = 1:1000
    for ii = 1:size(a,1)
        if strcmp(a(ii).name, strcat(prefix,num2str(jj),suffix))
            found = true;
            break
        else
            found = false;
        end 
    end
    
    if ~found
        loadj = jj-1;
        break
    end
end

clear j
% create bins
tmat = zeros(10,99);
% load file
for ii = 1:10
    load(strcat(prefix,num2str(ii),suffix),'rmyo','fc','nbead','ipt','ifor','rbead')
    if ii == 1
        circ_ring_noplot
        xbin = linspace(0,circring,100); 
    end
    tension_noplot
    for jj = 1:length(xbin)-1
        if any(and(xfit>xbin(jj), xfit<xbin(jj+1)))
            tmat(ii,jj) = mean(tens(and(xfit>xbin(jj), xfit<xbin(jj+1))));
        end
    end
    clear xfit
    clear tens
end