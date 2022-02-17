% change from original: also do sum intensity of Myo2 and Myp2
% change from v2: use gaussian instead of flat disk
% v3 and v4 skipped
% figure
% change from v5: add calculation of sep_myo_myp and bound_frac
% change from v6: calculate 80% interpercentile of actin 
figure
%% project components onto z-rho plane
[~,rho_myp] = cart2pol(rmyp(1,:),rmyp(2,:));
z_myp = rmyp(3,:);
[~,rho_myo] = cart2pol(rmyo(1,:),rmyo(2,:));
z_myo = rmyo(3,:);
[~,rho_act] = cart2pol(rbead(1,:),rbead(2,:));
z_act = rbead(3,:);
%% initialize the image plane
side_length = 0.4;
n_side = 100; % 100 squares, 101 X 101 lattice
im = zeros(n_side+1,n_side+1,3);
act_mat = false(n_side+1);
d_grid = side_length / n_side;
r_latt = linspace((r_ring+0.03-side_length),(r_ring+0.03),n_side+1);
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
    imone = imone + 2.^(-temp / rcapmyo_short^2);
end
imone = imone / max(imone(:));
im(:,:,1) = imone;
%% Myp2
imone = zeros(n_side+1,n_side+1);
for i = 1:size(rmyp,2)
    r = rho_myp(i);
    z = z_myp(i);
    temp = (R_latt - r).^2 + (Z_latt - z).^2;
    imone = imone + 2.^(-temp / rcapmyp^2);
end
if max(imone(:))~=0
    imone = imone / max(imone(:)) * .7;
end
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
% im = .8 * im;
im(:,:,1) = im(:,:,1) + act_mat*10000;
im(:,:,2) = im(:,:,2) + act_mat*10000;
im(:,:,3) = im(:,:,3) + act_mat*10000;
%% membrane
temp = r_latt-r_ring;
[~, temp] = min(abs(temp));
temp = n_side + 1 - temp;
im(temp,:,1) = 222/255;
im(temp,:,2) = 184/255;
im(temp,:,3) = 135/255;
imshow(im)

%% count what fraction of actin is bundled
    rho_myo = median(rho_myo);
    rho_myp = median(rho_myp);
    sep_myo_myp = rho_myo - rho_myp
    temp = [rbead(3,:);rho_act] - [median(rmyo(3,:));rho_myo];
    temp = sqrt(sum(temp.*temp));
    bound_frac = sum(temp < rcapmyo_short) / numel(temp);
    act_spread = prctile(z_act,95) - prctile(z_act,5)