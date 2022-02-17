% generate cylinder
r_cyl = r_ring;
[x,y,z] = cylinder;
surf(x*r_cyl, y*r_cyl, z*4-2, 'edgecolor','none','FaceColor',[255,153,51]/255,'FaceAlpha',.1)
hold on
% basic sphere
[x,y,z] = sphere;

% choose a particle
ipar = 1;
ipar2 =1;    

%% plot myp2 clusters as big cyan circles (not filled)
rad_myp = 0.13;
hold on
for i = 1:size(rmyp,2)
   surf(x*rad_myp+rmyp(1,i),y*rad_myp+rmyp(2,i),z*rad_myp+rmyp(3,i),'edgecolor','none','FaceColor','cyan','FaceAlpha',.2)
end
%% plot myo2 clusters as green circles
rad_myo = 0.066;
hold on
for i = 1:size(rmyo,2)
   surf(x*rad_myo+rmyo(1,i),y*rad_myo+rmyo(2,i),z*rad_myo+rmyo(3,i),'edgecolor','none','FaceColor','green','FaceAlpha',1)
end

%% plot first beads on a filament (formin) as a blue circle
rad_for = 0.004;
hold on
for i = 1:numel(ifor)
   surf(x*rad_for+rbead(1,ifor(i)),y*rad_for+rbead(2,ifor(i)),z*rad_for+rbead(3,ifor(i)),'edgecolor','none','FaceColor','blue','FaceAlpha',.5)
end  
    %% plot filament beads as small grey circles connected with thin lines
    % index of beads that are not formin dimers
    inotfor = 1:length(rbead);
    inotfor = inotfor(~ismember(inotfor,ifor));
        inotfor(inotfor == ipar) = [];
%     scatter3(rbead(1,inotfor),rbead(2,inotfor),rbead(3,inotfor),50,0.5 * [1 1 1],'fill')
    for i = 1:length(ifor)
        xx = rbead(1,ifor(i):ipt(i));
        yy = rbead(2,ifor(i):ipt(i));
        zz = rbead(3,ifor(i):ipt(i));
        plot3(xx, yy, zz,'Color',0.5*[1 1 1])
    end
    
     %% plot selected particle as a red dot
    rcat = [rbead, rmyo, rmyp];
%     scatter3(rcat(1,ipar),rcat(2,ipar),rcat(3,ipar),100,'red','fill')
%     scatter3(rcat(1,ipar2),rcat(2,ipar2),rcat(3,ipar2),100,'black','fill')
    %% plot crosslinkers as green lines
    % pairs of indices of crosslinkers
    [ix, jx] = find(xmat);
    for ii = 1:length(ix)
        i = ix(ii);
        j = jx(ii);
        plot3([rbead(1,i),rbead(1,j)],[rbead(2,i), rbead(2,j)],[rbead(3,i),rbead(3,j)],'Color',[0 1 0])
    end
    axis square;
    axis([-2,2,-2,2])
    axis equal
    view(0,90)
    hold off
    