% change from original: use weighted least squares
% load wtfc3_3_1sec
dir = '.';
prfx = 'gm5';
it = 120;
nval = 15;
i0 = 15;
run_max = 200;
ilist = i0:nval:run_max;

unit = 'sec';
load([dir '/' prfx '_1_00' unit '.mat']) 
%load(['wtfc3/wtfc3_8_10sec.mat']) % this one is good
% load([wtfc3/wtfc3_1_10sec.mat']) % this one is good


vmyo0 = 0.24;
bdmax = 60;
mean_tens = nan(length(ilist), bdmax);
sem_tens = nan(length(ilist), bdmax);
for idx= 1:length(ilist)
    irun = ilist(idx);
    if(exist([dir '/' prfx '_' num2str(irun) '_' num2str(it) unit '.mat'])==0)
        disp(['no file ' dir '/' prfx '_' num2str(irun) '_' num2str(it) unit '.mat'])
        continue
    end
    disp(['irun = ' num2str(irun)])
    load([dir '/' prfx '_' num2str(irun) '_' num2str(it) unit '.mat'])
    %% Calculate tension
    get_segtens           
    segtens_scalar1 = sqrt(sum(segtens.^2));
    rdiff = get_rdiff(rbead,1,ipt);
    temp = circshift(rdiff,1,2);
    segtens_scalar = -sum(temp.*segtens);

    beadind = zeros(1,length(segtens));
    for jseg = length(segtens):-1:1
        if any(ipt==jseg)
            beadind(jseg) = 1;
        else
            beadind(jseg) = beadind(jseg+1) + 1;
        end
    end
    % segtens_scalar_reverse = segtens_scalar(end:-1:1);
%     mean_tens(i, :) = nan(1,max(beadind));
%     sem_tens(i, :) = nan(1,max(beadind));
    for jseg = 1:bdmax
%         disp(jseg)
        n = sum(beadind==jseg);
        if n<10
%             disp('Not enough fils')
            continue
        end
        mean_tens(idx, jseg) = nanmean(segtens_scalar(beadind==jseg));
        sem_tens(idx, jseg) = nanstd(segtens_scalar(beadind==jseg))/sqrt(n-1);
    end
end

figure
box on
set(gca, 'Position',[0.2 0.2 0.7 0.7]);
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
hold on
xlabel('Distance from pointed end (\mum)')
ylabel('Filament tension (pN)')
errorbar((0:bdmax-1)/10, nanmean(mean_tens), nanmean(sem_tens), 'ko',...
    'LineWidth', 2, 'MarkerFaceColor', 'k')

%% theory plot
nhead = 16;
nx = length(rbead)*0.1/(2*pi*r_ring);
x1 = ((0:bdmax)/10)';
plot(x1, x1*nhead/nx*(1-(vpol+0.03)/vmyo0)*(length(rmyo)*fhead + length(rmyp)*fheadmyp)/...
    (2*pi*r_ring), 'Color', 'r', 'Linewidth', 2)
% axis([0, 1.5, 0, 35])

%% fit plot
% weighted least squares
mean_tens = nanmean(mean_tens);
sem_tens = nanmean(sem_tens);
if(any(isnan(sem_tens)))
    bdmax = find(isnan(sem_tens), 1)-2;
else
    bdmax = bdmax - 1;
end
x1 = ((0:bdmax)/10)';
% x1 = x1(1:ifil+1);
X = [ones(size(x1)),x1];
Y = mean_tens(1:bdmax+1)';
W = 1./sem_tens(1:bdmax+1);
a=lscov(X,Y,W);
plot(x1,X*a,'b-')
hold off

    saveas(gcf, [prfx '_tfil_i' num2str(i0) '.png']);
    saveas(gcf, [prfx '_tfil_i' num2str(i0) '.fig']);