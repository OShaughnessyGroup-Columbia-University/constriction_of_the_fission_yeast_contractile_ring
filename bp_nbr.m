function bpmat = bp_nbr (rbead, rmyp, ifor, ipt, rcapmyp)
if isempty(rbead)||isempty(rmyp)
    bpmat = logical.empty(size(rbead,2),size(rmyp,2));
    return
end
    

% %     % preallocate matrices that shows whether x and y positions of myosin
% %     % clusters and beads are within rcapmyp
% %     xtemp = zeros(length(rbead),size(rmyp,2));
% %     ytemp = xtemp;
    % find out possible neighbors whose x and y positions are within rcapmyp
    [xb, xm] = meshgrid(rbead(1,:), rmyp(1,:));
    [yb, ym] = meshgrid(rbead(2,:),rmyp(2,:));    
    [zb, zm] = meshgrid(rbead(3,:),rmyp(3,:));
    xd = xb - xm;
    yd = yb - ym;
    zd = zb - zm;
    rdsq = xd .* xd + yd .* yd + zd .* zd;
    temp = rdsq < rcapmyp * rcapmyp;
    bpmat = temp';
    
% %     
% %     
% %     for i = 1:length(rbead)
% %         xtemp(i,:) = abs(rbead(1,i)-rmyp(1,:)) < rcapmyp;
% %         ytemp(i,:) = abs(rbead(2,i)-rmyp(2,:)) < rcapmyp;
% %     end
% %     % make them sparse matrices
% %     xtemp = sparse(xtemp);
% %     ytemp = sparse(ytemp);
% %     % possible neighbors whose x and y coordinates are within rcapmyp
% %     bpmat = and(xtemp,ytemp);
       
    


% %     [beadind,myoind,value] = find(bpmat);
% %     % calculate distances for these possible neighbors and see if they are
% %     % within rcapmyp
% %     for i = 1:length(beadind)
% %         if norm(rbead(:,beadind(i)) - rmyp(:,myoind(i))) > rcapmyp 
% %             value(i)=false;
% %         end
% %     end
% %     bpmat = sparse(beadind,myoind,value,size(rbead,2),size(rmyp,2));
    
    
    
    % if more than one beads on the same filament are interacting with the same
    % myosin cluster, only the one that is closer to the pointed end but is
    % not the pointed itself counts (for simplicity)
    bpmat(ipt,:) = false;
    temp = ~circshift(bpmat,[-1,0]);
    bpmat = and(bpmat, temp);
    
% %     
% %     for j = 1:size(rmyp,2) 
% %         % boolean vector that contains all beads interacting with myosin #j
% %         bj = bpmat(:,j);
% %         % indices of all beads interacting with myosin #j
% %         beadindj = find(bj);
% %         for i = beadindj'
% %             if and(any(beadindj == i+1),... % if the next bead is also interacting with myosin #j
% %                      ~any(ipt == i)) % and this bead isn't a pointed end
% %                  if ~any(ipt == i+1) % if next bead isn't a pointed end
% %                     bpmat(i,j) = false;
% %                  else                % if next bead is a pointed end
% %                     bpmat(i+1,j) = false;
% %                  end
% %             elseif any(ipt == i)
% %                 bpmat(i,j) = false; % if this bead is a pointed end
% %             end
% %         end
% %     end

%     if select4      % if myp2 interacts with 4 actin filaments only
%         off_prob = dt / atp_cycle;  % probability that one myp2 head goes off the interacting actin filament in this time step
%         for i = 1:numel(bancm)  % loop over every myosin
%             if ~bancm(i)        % if this myosin is unanchored, i.e. myp2
%                 if sum(bpmat(:,i)) > 4  % if this myp2 is interacting with more than 4 actin filaments, need to select 4
%                     possible_ind = find(bpmat(:,i))';
%                     preferred_ind = find(bmmat_trial(:,i))';
%                     choice = [];
%                     for j = 1:numel(preferred_ind)
%                         if 1-rand < off_prob      % if this myp2 head will not go off
%                             temp = find(abs(possible_ind - preferred_ind(j)) < 2);
%                             choice = [choice, temp];
%                         end
%                     end
%                     if numel(choice) < 4                % if we haven't found 4 preferred filaments
%                         ind = 1:numel(possible_ind);
%                         ind(choice) = [];               % indices of those that are not already chosen
%                         temp = randperm(numel(ind), 4-numel(choice));
%                         choice = [choice, temp];
%                     end
%                     temp = possible_ind(choice); % indices of the 4 selected filaments to interact with this myp2
%                     bpmat(:,i) = false;
%                     bpmat(temp,i) = true;
%                 end
%             end
%         end
%     end
    if isempty(bpmat)
        bpmat = logical.empty(size(rbead,2),size(rmyp,2));
    else
        bpmat(ifor,:) = false;
    end
end