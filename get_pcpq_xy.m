function [pcpq, cvec, cdot] = get_pcpq_xy(rbead,rmyo,rmyp,ifor,fmmat,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo)
% will be run repeatedly
% constraints
% 1. x^2 + y^2 = r^2 for anchored components
% 2. fixed size between adjacent beads
% 3. formin and myosin in the same node fixed to the same theta and z

r_for_sq = (r_ring - d_for)^2;
r_myo_sq = (r_ring - d_myo)^2;

dbead_first_sq = dbead_first .* dbead_first;
    % concatenate x and y of rbead and rmyo
    xyzcat = [rbead,rmyo];
    % concatenate bancf and bancm
    banccat = [bancf, bancm];
    % number of anchored components
    nanc = sum(banccat);
    % number of actin beads (except formins)
    nact = size(rbead,2) - numel(ifor);
    % preallocate
    pcpq1 = zeros(nanc,3*(size(rbead,2)+size(rmyo,2)+size(rmyp,2)));
    c1 = zeros(1,nanc);
        cd1 = c1;
    c2 = zeros(1,nact);
        cd2 = c2;
    row2 = [];
    col2 = [];
    val2 = [];
% % %     row3 = [];
% % %     col3 = [];
% % %     val3 = [];
    %% anchored particles are fixed on a circle (plasma membrane)
    % indices of anchored particles
    ianc = 1:size(xyzcat,2);
    ianc = ianc(banccat);
    % form a 2 by nanc array containing 2x and 2y 
    % the factor 2 is unnecessary but kept for easier debugging
    twoxyz = 2 * xyzcat(:,ianc);
    % constraints for anchored particles
    for i = 1:nanc;
        pcpq1(i,3*(ianc(i)-1)+1) = twoxyz(1,i);
        pcpq1(i,3*(ianc(i)-1)+2) = twoxyz(2,i);
        pcpq1(i,3*(ianc(i)-1)+3) = 0*twoxyz(3,i);
        if i <= sum(bancf)
            c1(i) = xyzcat(1:2,ianc(i))' * xyzcat(1:2,ianc(i)) - r_for_sq;
        else
            c1(i) = xyzcat(1:2,ianc(i))' * xyzcat(1:2,ianc(i)) - r_myo_sq;
        end
        cd1(i) = 2 * r_ring * vs_abs;
    end
    
    %% distances between adjacent actin beads are fixed
    % indices of actin beads, excluding formin
    iact = 1:size(rbead,2);
    iact(ifor) = [];
    % form an array containing xy_this - xy_previous
    xyzdiff = rbead(:,:) - circshift(rbead(:,:),[0,1]);
    % exclude formin and multiply by 2
    twoxyzdiff = 2 * xyzdiff(:,iact);
    temp = sum(twoxyzdiff .* twoxyzdiff);
    c2 = .25 * temp - dbead_sq;
    % constraints for fixed actin subunit length
    for i = 1:nact        
        row2 = [row2; i;i;i;i;i;i];
        col2 = [col2; 3*(iact(i)-1)+1;3*(iact(i)-1)+2;3*iact(i);3*(iact(i)-2)+1;3*(iact(i)-2)+2;3*iact(i)-3]; 
        val2 = [val2; twoxyzdiff(1,i); twoxyzdiff(2,i);twoxyzdiff(3,i); -twoxyzdiff(1,i); -twoxyzdiff(2,i);-twoxyzdiff(3,i)];
        cd2(i) = 0;
    end
    
    % indices of first non-formin actin bead
    ind = ifor + 1;
    % change to boolean
    temp = false(1,size(rbead,2));
    temp(ind) = true;
    % exclude formins in the boolean vector
    temp(ifor) = [];
    % change to indices
    ind = find(temp);
    % change c2 and cd2 for the first non-formin actin beads
    for j = 1:numel(ind)
        i = ind(j);
        c2(i) = .25 * twoxyzdiff(:,i)' * twoxyzdiff(:,i) - dbead_first_sq(j);
        cd2(i) = - 2 * dbead_first(j) * vpol;
    end
    
    
    %% alpha-actinin crosslinks fix distance between crosslinked beads
% % %     pcpq3 = [];
% %     % if there is any crosslinker
% %     if sum(xmat(:)) > 0
% %         % take only the upper triangular part of xmat, to avoid double
% %         % counting
% %         xmatup = triu(xmat);
% %         % pairs of beads crosslinked by alpha-actinin
% %         [r,c] = find(xmatup);
% %         
% %         for i = 1:length(r)
% %             row3 = [row3; i;i;i;i;i;i];
% %             col3 = [col3; 3*(r(i)-1)+1; 3*(r(i)-1)+2; 3*r(i); 3*(c(i)-1)+1; 3*(c(i)-1)+2; 3*c(i)];
% %             val3 = [val3; rbead(:,r(i))-rbead(:,c(i)); rbead(:,c(i)) - rbead(:,r(i))];
% %             c3(i) = (rbead(:,r(i))-rbead(:,c(i)))' * (rbead(:,r(i))-rbead(:,c(i)));
% %         end
% %         
% %         wrong!! pcpq3 = sparse(row3,col3,val3,sum(xmatup(:)),3*(length(rbead)+length(rmyo)));
% %     else
% %         pcpq3 = [];
% %         c3 = [];
% %     end
    %% paired formin and myosin have the same theta and z
    % c3 = z_for - z_myo
    [i,j] = find(fmmat);
    c3 = nan(1,numel(i));
    cd3 = zeros(1,numel(i));
    pcpq3 = zeros(numel(i),3*(size(rbead,2)+size(rmyo,2)+size(rmyp,2)));
    for temp = 1:numel(i)
        ind_for = ifor(i(temp));
        z_for = rbead(3,ind_for);
        ind_myo  = j(temp);
        z_myo = rmyo(3,ind_myo);
        c3(temp)=z_for-z_myo;
        pcpq3(temp,3*ind_for) = 1;
        pcpq3(temp,3*(size(rbead,2)) + 3*ind_myo) = -1;
    end
    % c4 is set to reflect th_for == th_myo
    % in practice, c4 = x_for * y_myo - x_myo * y_for
    c4 = nan(1,numel(i));
    cd4 = zeros(1,numel(i));
    pcpq4 = zeros(numel(i),3*(size(rbead,2)+size(rmyo,2)+size(rmyp,2)));
    for temp = 1:numel(i)
        ind_for = ifor(i(temp));
        x_for = rbead(1,ind_for);
        y_for = rbead(2,ind_for);
        ind_myo  = j(temp);
        x_myo = rmyo(1,ind_myo);
        y_myo = rmyo(2,ind_myo);
        c4(temp)=  x_for * y_myo - x_myo * y_for;
        pcpq4(temp,3*ind_for-2) = y_myo;   % partial c/partial x_for = y_myo
        pcpq4(temp,3*ind_for-1) = -x_myo;  % partial c/partial y_for = - x_myo
        pcpq4(temp,3*(size(rbead,2)) + 3*ind_myo - 2) = - y_for;  % partial c/partial x_myo = - y_for
        pcpq4(temp,3*(size(rbead,2)) + 3*ind_myo - 1) = x_for;  % partial c/partial y_myo = x_for      
    end
    
    % concatenate the final matrix
    pcpq1 = sparse(pcpq1);
    pcpq2 = sparse(row2,col2,val2,nact,3*(size(rbead,2)+size(rmyo,2)+size(rmyp,2)));
    pcpq3 = sparse(pcpq3);
    pcpq4 = sparse(pcpq4);
    pcpq = [pcpq1;pcpq2;pcpq3;pcpq4];
    cvec = [c1, c2, c3, c4];
    cdot = [cd1, cd2, cd3, cd4];
end