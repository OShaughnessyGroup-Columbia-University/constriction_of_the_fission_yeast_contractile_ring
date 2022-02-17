% change from original: use weighted least squares
% load wtfc3_3_1min
% load(['wtfc3/wtfc3_11_2min.mat']) % this one is good
% load(['wtfc3/wtfc3_8_10min.mat']) % this one is good
load(['wtfc3/wtfc3_3_1min.mat']) % this one is good
% load(['wtfc3/wtfc3_5_10min.mat']) % this one is good
%% Calculate tension
get_segtens                    
segtens_scalar = sqrt(sum(segtens.^2));
beadind = zeros(1,length(segtens));
for i = length(segtens):-1:1
    if any(ipt==i)
        beadind(i) = 1;
    else
        beadind(i) = beadind(i+1) + 1;
    end
end
% segtens_scalar_reverse = segtens_scalar(end:-1:1);
mean_tens = nan(1,max(beadind));
sem_tens = nan(1,max(beadind));
for i = 1:max(beadind)
    mean_tens(i) = mean(segtens_scalar(beadind==i));
    n = sum(beadind==i);
    sem_tens(i) = std(segtens_scalar(beadind==i))/sqrt(n);
end
figure
% plot((1:max(beadind))/10,mean_tens,'ko')
hold on
errorbar((1:max(beadind))/10,mean_tens,sem_tens, 'ko', 'LineWidth', 2)
% weighted least squares
x1 = ((0:max(beadind))/10)';
x1 = x1(1:26);
X = [ones(size(x1)),x1];
Y = mean_tens(1:26)';
W = 1./sem_tens(1:26);
a=lscov(X,Y,W);
% plot(x1,X*a,'k-')
% NOTE: fudge factor needed. Is this compatible with predicted total ring
% tension?
plot(x1, 0.5*x1*(length(rmyo)*fhead + length(rmyp)*fheadmyp)/(2*pi*r_ring),...
    'Color', 'r', 'Linewidth', 2)
% plot(x1, 0.85*(1-28/240)*x1*(length(rmyo)*fhead + length(rmyp)*fheadmyp)/(2*pi*r_ring))
hold off
axis([0,1.5,0,45])
xlabel('Distance from the pointed end (\mum)')
ylabel('Filament tension (pN)')
box on
set(gca, 'linewidth', 2)
set(gca, 'fontsize', 30)
set(gca, 'position', [0.2 0.2 0.7 0.7])