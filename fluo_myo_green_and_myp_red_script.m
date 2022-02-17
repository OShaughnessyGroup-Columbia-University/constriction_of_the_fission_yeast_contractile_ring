%% myp2 part

pixel_size = 0.01;
sigma = 0.600*.61/1.4/2; % width of the Gaussian point-spread function. NA = 1.4 for Laplante 3 myosin 
inv_s_sq = (sigma+0.1)^-2;
rplot = rmyp(:,abs(rmyp(3,:))<0.3);

theta = linspace(0,2*pi,9);
    theta(end) = [];

left = -2;
right = 2;
top = 2;
bottom = -2;

x_vec = left:pixel_size:right;
y_vec = top:-pixel_size:bottom;
fluo = zeros(numel(y_vec), numel(x_vec));
[x,y] = meshgrid(x_vec,y_vec);

fluo = zeros(numel(y_vec), numel(x_vec));
for i = 1:size(rplot,2)
    this_rmyo = rplot(:,i);
    dx = x - this_rmyo(1);
    dy = y - this_rmyo(2);
    dr_sq = dx .* dx + dy .* dy;
    r2 = this_rmyo .* this_rmyo;
    r2 = r2(1)+r2(2);
    if r2 > 1
    fluo = fluo + exp(-dr_sq .* .5 * inv_s_sq);
    else        
    fluo = fluo + exp(-dr_sq .* .5 * inv_s_sq);
    end
end

im = zeros(size(fluo,1),size(fluo,2),3);
im(:,:,1) = fluo /max(fluo(:))640;

% imshow(imrotate(im,-90))

%% myo2 part

inv_s_sq = (sigma+0.05)^-2;
rplot = rmyo(:,abs(rmyo(3,:))<0.3);

theta = linspace(0,2*pi,9);
    theta(end) = [];

fluo = zeros(numel(y_vec), numel(x_vec));
for i = 1:size(rplot,2)
    this_rmyo = rplot(:,i);
    dx = x - this_rmyo(1);
    dy = y - this_rmyo(2);
    dr_sq = dx .* dx + dy .* dy;
    r2 = this_rmyo .* this_rmyo;
    r2 = r2(1)+r2(2);
    if r2 > 1
    fluo = fluo + exp(-dr_sq .* .5 * inv_s_sq);
    else
    fluo = fluo + exp(-dr_sq .* .5 * inv_s_sq);
    end
end

im(:,:,2) = fluo /max(fluo(:));

imshow(imrotate(im,0))
% 
% %% plot radial distribution of myo2 and myp2
% figure 
% 
% rho = sqrt(x.*x + y.*y);
% nbins = 100;
% immyo = im(:,:,1);
% immyp = im(:,:,2);
% [histw_myo, intervals_myo] = histwc(rho(:), immyo(:), nbins);
% histw_myo = histw_myo / max(histw_myo);
% plot(intervals_myo, histw_myo, '-r','LineWidth',3)
% hold on
% [histw_myp, intervals_myp] = histwc(rho(:), immyp(:), nbins);
% histw_myp = histw_myp / max(histw_myp);
% plot(intervals_myp, histw_myp, '-g','LineWidth',3)
% hold off 
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16)
% axis([0,1.5,0,1])
% xlabel('Radial position (\mum)','FontSize',20)
% ylabel('Normalized fluorescence intensity','FontSize',20)
% legend('Myo2p','Myp2p','Location','northwest')
% 
% [~,radmyo2]=cart2pol(rmyo(1,:),rmyo(2,:));
% [~,radmyp2]=cart2pol(rmyp(1,:),rmyp(2,:));
% mean(radmyo2)-mean(radmyp2)