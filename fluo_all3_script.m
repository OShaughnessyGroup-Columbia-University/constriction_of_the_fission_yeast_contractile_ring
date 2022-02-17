%% myp2 part

pixel_size = 0.01;
sigma = 0.09;
inv_s_sq = (sigma+0.1/3)^-2;
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
im(:,:,2) = fluo /max(fluo(:));

% imshow(imrotate(im,-90))

%% myo2 part

inv_s_sq = (sigma+0.05/3)^-2;
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

im(:,:,1) = fluo /max(fluo(:));



%% actin part
inv_s_sq = sigma^-2;
rplot = rbead(:,abs(rbead(3,:))<0.3);

theta = linspace(0,2*pi,9);
    theta(end) = [];

left = -2;
right = 2;
top = 2;
bottom = -2;

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

im(:,:,3) = fluo /max(fluo(:));

imshow(imrotate(im,0))