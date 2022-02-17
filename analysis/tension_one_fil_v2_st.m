% change from original: use weighted least squares
% load wtfc3_3_1min
load(['wtfc3/wtfc3_11_00min.mat']) % this one is good
load(['wtfc3/wtfc3_11_10min.mat']) % this one is good
vmyo0=0.24;
%% Calculate tension
get_segtens           
segtens_scalar1 = sqrt(sum(segtens.^2));

rdiff = get_rdiff(rbead,1,ipt);
temp = circshift(rdiff,1,2);
segtens_scalar = -sum(temp.*segtens);

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
    n = sum(beadind==i);
    if n<60
        continue
    end
    mean_tens(i) = mean(segtens_scalar(beadind==i));
    sem_tens(i) = std(segtens_scalar(beadind==i))/sqrt(n);
    sd_tens(i) = std(segtens_scalar(beadind==i));
end
ifil = round(mean(beadind));
figure
% plot((1:max(beadind))/10,mean_tens,'ko')
hold on
errorbar((1:ifil)/10, mean_tens(1:ifil), sem_tens(1:ifil), 'ko',...
    'LineWidth', 2, 'MarkerFaceColor', 'k')
% weighted least squares
x1 = ((0:15)/10)';
% x1 = x1(1:ifil+1);
% X = [ones(size(x1)),x1];
% Y = mean_tens(1:25)';
% W = 1./sem_tens(1:25);
% a=lscov(X,Y,W);
% plot(x1,X*a,'k-')

nhead = 16;
box on
set(gca, 'Position',[0.2 0.2 0.7 0.7]);
set(gca, 'Linewidth', 2)
set(gca, 'FontSize', 30)
ft = fittype('a*x');
[lin gof] = fit((0.1:0.1:1.4)', mean_tens(1:14)', ft)
nx = length(rbead)*0.1/(2*pi*r_ring);
plot(x1, nhead/nx*(1-vpol/vmyo0)*x1*(length(rmyo)*fhead + length(rmyp)*fheadmyp)...
    /(2*pi*r_ring), 'Color', 'b', 'Linewidth', 2)
plot(x1, lin(x1), 'Color', 'r', 'Linewidth', 2)
hold off
axis([0, 1.5, 0, 35])
xlabel('Distance from pointed end (\mum)')
ylabel('Filament tension (pN)')

%% chi squared
% res = mean_tens(1:ifil) - x1(2:1+ifil)'*nhead/nx*(1-vpol/vmyo0)*...
%     (length(rmyo)*fhead + length(rmyp)*fheadmyp)/(2*pi*r_ring)
% res2 = res.^2;% - sem_tens(1:ifil).^2;
% chi2 = sum(res2./sd_tens(1:ifil).^2)
res = mean_tens(1:ifil) - lin(x1(2:15))';
res2 = res.^2;% - sem_tens(1:ifil).^2;
chi2 = sum(res2./sd_tens(1:ifil).^2)


% res = mean_tens(1:ifil) - x1(2:1+ifil)'*nhead/nx*(1-vpol/vmyo0)*...
%     (length(rmyo)*fhead + length(rmyp)*fheadmyp)/(2*pi*r_ring)
% res2 = res.^2;
% chi2 = sum(res2./sd_tens(1:ifil).^2)

