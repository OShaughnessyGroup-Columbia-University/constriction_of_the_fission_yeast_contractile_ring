
% choose a particle
ipar = 1;
ipar2 =2;    
    % change bf (all beads) to bf (formin only)
    bf = bancf(ifor);
    %% plot myo2 clusters as big orange circles (filled - anchored)
    imyo = 1:size(rmyo,2);
    % plot unanchored
    iunanc = imyo(~bancm);
        iunanc(iunanc == ipar - size(rbead,2)) = [];
    scatter3(rmyo(1,iunanc),rmyo(2,iunanc),rmyo(3,iunanc),100,'red')
    hold on
    % plot anchored
    ianc = imyo(bancm);
        ianc(ianc == ipar - size(rbead,2)) = [];
    scatter3(rmyo(1,ianc),rmyo(2,ianc),rmyo(3,ianc),100,'red','fill')
    %% plot myp2 clusters as big orange circles (not filled) connected with orange lines
    % might be confused with unanchored myo2. However currently all myo2
    % are anchored
    imyp = 1:size(rmyp,2);
    % plot unanchored
    scatter3(rmyp(1,:),rmyp(2,:),rmyp(3,:),100,'green')
%     for i = 1:(size(rmyp,2)/2)
%         xx = rmyp(1,2*i-1:2*i);
%         yy = rmyp(2,2*i-1:2*i);
%         zz = rmyp(3,2*i-1:2*i);
%         plot3(xx, yy, zz,'Color',[255,165,0]/255)
%     end
    %% plot first beads on a filament (formin) as a blue circle (filled - anchored)
    % next bead of ipt is the formin bead
    ifor = ipt + 1;
    % delete the last bead (which does not exist) and add the first bead
    ifor = ifor(1:end-1);
    ifor = [1,ifor];
    % plot unanchored
    iunanc = ifor(~bf);
        iunanc(iunanc == ipar) = [];
    scatter3(rbead(1,iunanc),rbead(2,iunanc),rbead(3,iunanc),50,'blue')
    % plot anchored
    ianc = ifor(bf);
        ianc(ianc == ipar) = [];
    scatter3(rbead(1,ianc),rbead(2,ianc),rbead(3,ianc),5,'blue','fill')
    %% plot filament beads as small grey circles connected with thin lines
    % index of beads that are not formin dimers
    inotfor = 1:size(rbead,2);
    inotfor = inotfor(~ismember(inotfor,ifor));
        inotfor(inotfor == ipar) = [];
    scatter3(rbead(1,inotfor),rbead(2,inotfor),rbead(3,inotfor),50,0.5 * [1 1 1],'fill')
    for i = 1:length(ifor)
        xx = rbead(1,ifor(i):ipt(i));
        yy = rbead(2,ifor(i):ipt(i));
        zz = rbead(3,ifor(i):ipt(i));
        plot3(xx, yy, zz,'Color',0.5*[1 1 1])
    end
    
     %% plot selected particle as a red dot
    rcat = [rbead, rmyo, rmyp];
    scatter3(rcat(1,ipar),rcat(2,ipar),rcat(3,ipar),100,'magenta','fill')
    scatter3(rcat(1,ipar2),rcat(2,ipar2),rcat(3,ipar2),100,'black','fill')
    %% plot crosslinkers as green lines
    % pairs of indices of crosslinkers
    [ix, jx] = find(xmat);
    for ii = 1:length(ix)
        i = ix(ii);
        j = jx(ii);
        plot3([rbead(1,i),rbead(1,j)],[rbead(2,i), rbead(2,j)],[rbead(3,i),rbead(3,j)],'Color',[0 1 0])
    end
    %% plot the cell wall
    radius = r_ring; 
    centerX = 0;
    centerY = 0;
    rectangle('Position',[centerX - radius, centerY - radius, radius*2, radius*2],...
        'Curvature',[1,1]);
    axis square;
   
    axis equal
    view(0,90)
    