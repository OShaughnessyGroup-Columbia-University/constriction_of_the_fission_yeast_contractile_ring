% change from original: generate an additional plot of typical component
% organization

%% typical component organization
hold on
% basic sphere
[x,y,z] = sphere;
% Myp2p
rad_myp = rcapmyp;
[~,rho_myp] = cart2pol(rmyp(1,:),rmyp(2,:));
rho_myp = median(rho_myp);
circles(median(rmyp(3,:)),rho_myp,rad_myp,'facecolor','green','edgecolor','green','facealpha',.3)
% circles(median(rmyp(3,:)),rho_myp,rad_myp/2,'facecolor','none','edgecolor','green')
% Myo2p
rad_myo = rcapmyo_short;
[~,rho_myo] = cart2pol(rmyo(1,:),rmyo(2,:));
rho_myo = median(rho_myo);
circles(median(rmyo(3,:)),rho_myo,rad_myo,'facecolor','red','edgecolor','red','facealpha',.3)
% circles(median(rmyo(3,:)),rho_myo,rad_myo/2,'facecolor','none','edgecolor','red')
    %% plot filament beads as small grey circles connected with thin lines
    rad_act = 0.001;
    [~,rho_act] = cart2pol(rbead(1,:),rbead(2,:));
    randind = randi(size(rbead,2),1,2000);
    scatter(rbead(3,randind),rho_act(randind),10,'k.')
%     for i = randi(size(rbead,2),1,size(rbead,2))
%        surf(x*rad_act+rbead(3,i),y*rad_act+rho_act(i),z*rad_act+20,'edgecolor','none','FaceColor','black')
%     end
    % plot ceiling
    plot([-.3,.3],[r_ring,r_ring])
    % plot levels of Myo2p and Myp2p
%     plot([-.3,.3],[rho_myp,rho_myp],'g-')
%     plot([-.3,.3],[rho_myo,rho_myo],'r-')
    axis equal
    axis([-0.15,0.15,r_ring - 0.29,r_ring + 0.01])
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    axis off
    hold off
    sep_myo_myp = rho_myo - rho_myp;
    %% count what fraction of actin is bundled
    temp = [rbead(3,:);rho_act] - [median(rmyo(3,:));rho_myo];
    temp = sqrt(sum(temp.*temp));
    bound_frac = sum(temp < rcapmyo_short) / numel(temp);