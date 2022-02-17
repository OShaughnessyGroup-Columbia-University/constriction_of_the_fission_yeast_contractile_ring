% change from original: also do sum intensity of Myo2 and Myp2
% change from v2: use gaussian instead of flat disk
% v3 and v4 skipped
% figure
% change from v5: add calculation of sep_myo_myp and bound_frac
% change from v6: calculate 80% interpercentile of actin 
% change from v7: calculate 80% interpercentile of actin 
% load('delmyp7_2_10min.mat');
figure
%% initialize the image plane
side_length = 2*r_ring;
n_side = 400; % 100 squares, 101 X 101 lattice
im = zeros(n_side+1,n_side+1,3);
act_mat = false(n_side+1);
d_grid = side_length / n_side;
x_latt = linspace(-r_ring-0.03, r_ring+0.03, n_side+1);
y_latt = linspace(-r_ring-0.03, r_ring+0.03, n_side+1);
[X_latt, Y_latt] = meshgrid(x_latt,y_latt);
%Y_latt = Y_latt;
%% Myo2
imone = zeros(n_side+1,n_side+1);
for i = 1:size(rmyo,2)
    x = rmyo(1, i);
    y = rmyo(2, i);
    temp = (X_latt - x).^2 + (Y_latt - y).^2;
    imone = imone + 2.^(-temp / rcapmyo_short^2);
end
imone = imone / max(imone(:));
im(:,:,1) = imone;
% imshow(im)
%% Myp2
if(~isempty(rmyp))
    imone = zeros(n_side+1,n_side+1);
    for i = 1:size(rmyp,2)
        x = rmyp(1, i);
        y = rmyp(2, i);
        temp = (X_latt - x).^2 + (Y_latt - y).^2;
        imone = imone + 2.^(-temp / rcapmyp^2);
    end
    if max(imone(:))~=0
        imone = imone / max(imone(:)) * .7;
    end
    im(:,:,2) = im(:,:,2) + imone;
end
% % imshow(im)
%% actin
for i = 1:size(rbead,2)
    % only plot 10 percent
    if rand > .5
        continue
    end
    x = rbead(1, i);
    y = rbead(2, i);
    temp = (X_latt - x).^2 + (Y_latt - y).^2;
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
temp = X_latt.^2 + Y_latt.^2 - r_ring^2;
[mrow, mcol] = find(abs(temp) < 0.01);
% temp = n_side + 1 - temp;
for i=1:length(mrow)
    im(mrow(i),mcol(i),1) = 222/255;
    im(mrow(i),mcol(i),2) = 184/255;
    im(mrow(i),mcol(i),3) = 135/255;
end

%% Scale bar
lbar = .200;
nlbar = floor(lbar/side_length*n_side);
wbar = .020;
nwbar = floor(wbar/side_length*n_side);
im(50:50+nwbar, 350:350+nlbar, 1) = 10000;
im(50:50+nwbar, 350:350+nlbar, 2) = 10000;
im(50:50+nwbar, 350:350+nlbar, 3) = 10000;


%% show
imshow(im)