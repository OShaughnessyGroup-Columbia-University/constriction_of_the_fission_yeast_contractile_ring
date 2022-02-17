 % choose a particle
ipar = 1;
    
    % change bf (all beads) to bf (formin only)
    bf = bancf(ifor);
    %% plot myosin clusters as big orange circles (filled - anchored)
    imyo = 1:length(rmyo);
    % plot unanchored
    iunanc = imyo(~bancm);
        iunanc(iunanc == ipar - nbead) = [];
    scatter3(rmyo(1,iunanc),rmyo(2,iunanc),rmyo(3,iunanc),100,[255,165,0]/255)
    hold on
    % plot anchored
    ianc = imyo(bancm);
        ianc(ianc == ipar - nbead) = [];
    scatter3(rmyo(1,ianc),rmyo(2,ianc),rmyo(3,ianc),100,[255,165,0]/255,'fill')
    
    
     %% plot selected particle as a red dot
    rcat = [rbead, rmyo];
    scatter3(rcat(1,ipar),rcat(2,ipar),rcat(3,ipar),100,'red','fill')
    
    %% plot crosslinkers as green lines
    % pairs of indices of crosslinkers
    [ix, jx] = find(xmat);
    for ii = 1:length(ix)
        i = ix(ii);
        j = jx(ii);
        plot3([rbead(1,i),rbead(1,j)],[rbead(2,i), rbead(2,j)],[rbead(3,i),rbead(3,j)],'Color',[0 1 0])
    end
    
   
    axis equal
    view(0,90)