close all

% basic sphere
[x,y,z] = sphere;

% choose a particle
ipar = 1;
ipar2 =1;    

%% plot myp2 clusters as big cyan circles (not filled)
rad_myp = rcapmyp;
[~,rho_myp] = cart2pol(rmyp(1,:),rmyp(2,:));
hold on
for i = 1:size(rmyp,2)
   surf(x*rad_myp+rmyp(3,i),y*rad_myp+rho_myp(i),z*rad_myp,'edgecolor','none','FaceColor','green','FaceAlpha',10/size(rmyp,2))
end
%% plot myo2 clusters as green circles
rad_myo = 0.054;
[~,rho_myo] = cart2pol(rmyo(1,:),rmyo(2,:));
hold on
for i = 1:size(rmyo,2)
   surf(x*rad_myo+rmyo(3,i),y*rad_myo+rho_myo(i),z*rad_myo+10,'edgecolor','none','FaceColor','red','FaceAlpha',5/size(rmyo,2))
end

    %% plot filament beads as small grey circles connected with thin lines
    rad_act = 0.001;
    [~,rho_act] = cart2pol(rbead(1,:),rbead(2,:));
    for i = randi(size(rbead,2),1,size(rbead,2))
       surf(x*rad_act+rbead(3,i),y*rad_act+rho_act(i),z*rad_act+20,'edgecolor','none','FaceColor','black')
    end
    % plot ceiling
    plot([-.3,.3],[r_ring,r_ring])
    axis equal
    hold off
    
%% plot radial distributions
figure 

x_vec = left:pixel_size:right;
y_vec = bottom:pixel_size:top;
[x,y] = meshgrid(x_vec,y_vec);

rho = sqrt(x.*x + y.*y);
nbins = 100;
immyo = im(:,:,1);
immyp = im(:,:,2);
imact = im(:,:,3);
[histw_myo, intervals_myo] = histwc(rho(:), immyo(:), nbins);
histw_myo = histw_myo / max(histw_myo);
plot(intervals_myo, histw_myo, '-r','LineWidth',3)
hold on
[histw_myp, intervals_myp] = histwc(rho(:), immyp(:), nbins);
histw_myp = histw_myp / max(histw_myp);
plot(intervals_myp, histw_myp, '-g','LineWidth',3)
[histw_act, intervals_act] = histwc(rho(:), imact(:), nbins);
histw_act = histw_act / max(histw_act);
plot(intervals_act, histw_act, '-b','LineWidth',3)
hold off 
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
axis([0,1.5,0,1])
xlabel('Radial position (\mum)','FontSize',18)
ylabel('Normalized fluorescence intensity','FontSize',18)
legend('Myo2p','Myp2p','Actin','Location','northwest')

[~,radmyo2]=cart2pol(rmyo(1,:),rmyo(2,:));
[~,radmyp2]=cart2pol(rmyp(1,:),rmyp(2,:));
mean(radmyo2)-mean(radmyp2)