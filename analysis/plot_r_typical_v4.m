% change from original: also do sum intensity of Myo2 and Myp2
% change from v2: use white background naively
% v3 skipped
% figure
%% project components onto z-rho plane
[~,rho_myp] = cart2pol(rmyp(1,:),rmyp(2,:));
z_myp = rmyp(3,:);
[~,rho_myo] = cart2pol(rmyo(1,:),rmyo(2,:));
z_myo = rmyo(3,:);
[~,rho_act] = cart2pol(rbead(1,:),rbead(2,:));
z_act = rbead(3,:);
%% initialize the image plane
side_length = 0.3;
n_side = 100; % 100 squares, 101 X 101 lattice
im = zeros(n_side+1,n_side+1,3);
act_mat = false(n_side+1);
d_grid = side_length / n_side;
r_latt = linspace((r_ring+0.01-side_length),(r_ring+0.01),n_side+1);
z_latt = linspace(-side_length/2, side_length/2,n_side+1);
[R_latt, Z_latt] = meshgrid(r_latt,z_latt);
R_latt = flipud(R_latt');
Z_latt = Z_latt';
%% Myo2
imone = zeros(n_side+1,n_side+1);
for i = 1:size(rmyo,2)
    r = rho_myo(i);
    z = z_myo(i);
    temp = (R_latt - r).^2 + (Z_latt - z).^2;
    imone = imone + (temp < rcapmyo_short^2);
end
imone = imone / max(imone(:));
im(:,:,1) = imone;
%% Myp2
imone = zeros(n_side+1,n_side+1);
for i = 1:size(rmyp,2)
    r = rho_myp(i);
    z = z_myp(i);
    temp = (R_latt - r).^2 + (Z_latt - z).^2;
    imone = imone + (temp < rcapmyp^2);
end
imone = imone / max(imone(:));
im(:,:,2) = im(:,:,2) + imone;
%% actin
for i = 1:size(rbead,2)
    % only plot 10 percent
    if rand > .1
        continue
    end
    r = rho_act(i);
    z = z_act(i);
    temp = (R_latt - r).^2 + (Z_latt - z).^2;
    [row, col] = find(temp == min(temp(:)));
    act_mat(row,col) = true;
end
act_mat(1,:) = false;
act_mat(:,1) = false;
act_mat(end,:) = false;
act_mat(:,end) = false;
%% make background white
temp = sum(im,3);
[d1vec,d2vec] = find (temp < 0.01);
for i = 1:length(d1)
    d1 = d1vec(i);
    d2 = d2vec(i);
    im(d1,d2,1) = 1;
    im(d1,d2,2) = 1;
    im(d1,d2,3) = 1;
end
imshow(im)