%% %% change from v1: myosin does not interact with all actin beads within its capture radius. Instead, it selects 16 randomly to interact with.

function bpmat = bp_nbr_v2 (rbead, rmyp, ifor, ipt, rcapmyp)
if isempty(rbead)||isempty(rmyp)
    bpmat = logical.empty(size(rbead,2),size(rmyp,2));
    return
end
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
    
    % if more than one beads on the same filament are interacting with the same
    % myosin cluster, only the one that is closer to the pointed end but is
    % not the pointed itself counts (for simplicity)
    bpmat(ipt,:) = false;
    temp = ~circshift(bpmat,[-1,0]);
    bpmat = and(bpmat, temp);
    
    %% each myo2 node interacts with 16 actin filaments only
for i = 1:size(bpmat,2)  % loop over every myosin
    if sum(bpmat(:,i)) > 16  % if this myo2 is interacting with more than 16 actin filaments, need to select 4
        possible_ind = find(bpmat(:,i))';
        bpmat(:,i) = false;
        bpmat(randsample(possible_ind,16),i) = true;
    end
end
    %%
    if isempty(bpmat)
        bpmat = logical.empty(size(rbead,2),size(rmyp,2));
    else
        bpmat(ifor,:) = false;
    end
end