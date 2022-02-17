function bmmatpt = ptcap(rbead, rmyo, ipt, bancm, rcapmyo_u, rcapmyo_long, rcapmyo_short, bmmat)
% bmmatpt is the same as bmmat except at the pointed end
[row,col] = find(bmmat);
for i = 1:numel(row)
    thisrow = row(i);
    if any(ipt == i+1)          
        thiscol = col(i);
        dr = rbead(:,thisrow) - rmyo(:,thiscol);
        if ~bancm(thiscol)&&sum(dr.*dr) < rcapmyo_u*rcapmyo_u
            row(i) = thisrow + 1;
        elseif bancm(thiscol)
            xm = rmyo(1,thiscol);
            ym = rmyo(2,thiscol);
            zm = rmyo(3,thiscol);
            [thm,~] = cart2pol(xm,ym);
            r_hat = [cos(thm);sin(thm)];
            th_hat = [-sin(thm);cos(thm)];
            xd = rbead(1,thisrow) - xm;
            yd = rbead(2,thisrow) - ym;
            zd = rbead(3,thisrow) - zm;
            rd = [xd', yd'] * r_hat;
            thd = [xd', yd'] * th_hat;
            temp = (rd/rcapmyo_short).^2 + (thd/rcapmyo_long).^2 + (zd'/rcapmyo_short).^2;
            if temp < 1
                row(i) = thisrow + 1;
            end
        end
    end
end
bmmatpt = sparse(row,col,true(numel(row),1),size(rbead,2),size(rmyo,2));
end