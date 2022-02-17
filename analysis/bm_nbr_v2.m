%% change from v1: myosin does not interact with all actin beads within its capture radius. Instead, it selects 16 randomly to interact with.

function bmmat = bm_nbr_v2 (rbead, rmyo, ifor, ipt, rcapmyo_long, rcapmyo_short, fmmat)
if isempty(rbead)||isempty(rmyo)
    bmmat = logical.empty(size(rbead,2),size(rmyo,2));
    return
end
    

% %     % preallocate matrices that shows whether x and y positions of myosin
% %     % clusters and beads are within rcap
% %     xtemp = zeros(length(rbead),size(rmyo,2));
% %     ytemp = xtemp;
    % find out possible neighbors whose x and y positions are within rcap
    [xb, xm] = meshgrid(rbead(1,:), rmyo(1,:));
    [yb, ym] = meshgrid(rbead(2,:),rmyo(2,:));    
    [zb, zm] = meshgrid(rbead(3,:),rmyo(3,:));
    xd = xb - xm;
    yd = yb - ym;
    zd = zb - zm;
    rdsq = xd .* xd + yd .* yd + zd .* zd;
    temp = rdsq < rcapmyo_long * rcapmyo_long;
    bmmat = temp';
    
    for i = 1:size(rmyo,2)
        xm = rmyo(1,i);
        ym = rmyo(2,i);
        zm = rmyo(3,i);
        [thm,~] = cart2pol(xm,ym);
        r_hat = [cos(thm);sin(thm)];
        th_hat = [-sin(thm);cos(thm)];
        xd = rbead(1,:) - xm;
        yd = rbead(2,:) - ym;
        zd = rbead(3,:) - zm;
        rd = [xd', yd'] * r_hat;
        thd = [xd', yd'] * th_hat;
        temp = (rd/rcapmyo_short).^2 + (thd/rcapmyo_long).^2 + (zd'/rcapmyo_short).^2;
        bmmat(:,i) = (temp < 1);
    end
% %     
% %     
% %     for i = 1:length(rbead)
% %         xtemp(i,:) = abs(rbead(1,i)-rmyo(1,:)) < rcap;
% %         ytemp(i,:) = abs(rbead(2,i)-rmyo(2,:)) < rcap;
% %     end
% %     % make them sparse matrices
% %     xtemp = sparse(xtemp);
% %     ytemp = sparse(ytemp);
% %     % possible neighbors whose x and y coordinates are within rcap
% %     bmmat = and(xtemp,ytemp);
       
    


% %     [beadind,myoind,value] = find(bmmat);
% %     % calculate distances for these possible neighbors and see if they are
% %     % within rcap
% %     for i = 1:length(beadind)
% %         if norm(rbead(:,beadind(i)) - rmyo(:,myoind(i))) > rcap 
% %             value(i)=false;
% %         end
% %     end
% %     bmmat = sparse(beadind,myoind,value,size(rbead,2),size(rmyo,2));
    
    
    
    % if more than one beads on the same filament are interacting with the same
    % myosin cluster, only the one that is closer to the pointed end but is
    % not the pointed itself counts (for simplicity)
    bmmat(ipt,:) = false;
    temp = ~circshift(bmmat,[-1,0]);
    bmmat = and(bmmat, temp);
    
% %     
% %     for j = 1:size(rmyo,2) 
% %         % boolean vector that contains all beads interacting with myosin #j
% %         bj = bmmat(:,j);
% %         % indices of all beads interacting with myosin #j
% %         beadindj = find(bj);
% %         for i = beadindj'
% %             if and(any(beadindj == i+1),... % if the next bead is also interacting with myosin #j
% %                      ~any(ipt == i)) % and this bead isn't a pointed end
% %                  if ~any(ipt == i+1) % if next bead isn't a pointed end
% %                     bmmat(i,j) = false;
% %                  else                % if next bead is a pointed end
% %                     bmmat(i+1,j) = false;
% %                  end
% %             elseif any(ipt == i)
% %                 bmmat(i,j) = false; % if this bead is a pointed end
% %             end
% %         end
% %     end
%% each myo2 node interacts with 16 actin filaments only
for i = 1:size(bmmat,2)  % loop over every myosin
    if sum(bmmat(:,i)) > 16  % if this myo2 is interacting with more than 16 actin filaments, need to select 4
        possible_ind = find(bmmat(:,i))';
        bmmat(:,i) = false;
        bmmat(randsample(possible_ind,16),i) = true;
    end
end
    %% myosin does not bind or pull actin from its own node
    [i,j] = find(fmmat);
    for temp = 1:length(i)
        this_i = i(temp);
        this_j = j(temp);
        this_fil = ifor(this_i):ipt(this_i);
        bmmat(this_fil,this_j) = false;
    end
    
    %% if bmmat is empty
    if isempty(bmmat)
        bmmat = logical.empty(size(rbead,2),size(rmyo,2));
    end
end