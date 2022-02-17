function bpmatpt = ptcap_myp(rbead, rmyp, ipt, rcapmyp, bpmat)
% bmmatpt is the same as bmmat except at the pointed end
[row,col] = find(bpmat);
for i = 1:numel(row)
    thisrow = row(i);
    if any(ipt == i+1)          
        thiscol = col(i);
        dr = rbead(:,thisrow) - rmyp(:,thiscol);
        if sum(dr.*dr) < rcapmyp
            row(i) = thisrow + 1;
        end
    end
end
bpmatpt = sparse(row,col,true(numel(row),1),size(rbead,2),size(rmyp,2));
end