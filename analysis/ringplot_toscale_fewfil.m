 % choose a particle
ipar = 1;
ipar2 =1;    
    % change bf (all beads) to bf (formin only)
    bf = bancf(ifor);
    %% plot myp2 clusters as big orange circles (not filled) connected with orange lines
    % might be confused with unanchored myo2. However currently all myo2
    % are anchored
    imyp = 1:size(rmyp,2);
    % plot unanchored  
    if ~isempty(rmyp)
        circles(rmyp(1,:),rmyp(2,:),.1,'facecolor','none','edgecolor','green')
    end
    %% plot myo2 clusters as green circles (filled - anchored)
     imyo = 1:size(rmyo,2);
%     % plot unanchored
%     iunanc = imyo(~bancm);
%         iunanc(iunanc == ipar - size(rbead,2)) = [];
    %scatter3(rmyo(1,iunanc),rmyo(2,iunanc),rmyo(3,iunanc),100,[255,165,0]/255)
    hold on
    % plot anchored
    ianc = imyo(bancm);
        ianc(ianc == ipar - size(rbead,2)) = [];
    circles(rmyo(1,ianc),rmyo(2,ianc),.051,'facecolor','none','edgecolor','red')
%     scatter3(rmyo(1,ianc),rmyo(2,ianc),rmyo(3,ianc),100,[0,204,102]/255,'fill')

    %% plot first beads on a filament (formin) as a blue circle (filled - anchored)
    circles(rbead(1,ifor),rbead(2,ifor),.0035,'facecolor',[0,0,255]/255,'edgecolor','none')
    %% plot filament beads as small grey circles connected with thin lines
    % index of beads that are not formin dimers
    inotfor = 1:length(rbead);
    inotfor = inotfor(~ismember(inotfor,ifor));
        inotfor(inotfor == ipar) = [];
%     scatter3(rbead(1,inotfor),rbead(2,inotfor),rbead(3,inotfor),50,0.5 * [1 1 1],'fill')
    for i = randsample(length(ifor),5)'
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
        plot3([rbead(1,i),rbead(1,j)],[rbead(2,i), rbead(2,j)],[rbead(3,i),rbead(3,j)],'Color','cyan')
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
    