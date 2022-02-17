function ncross = crossDetect(pfx, irun, it)
% algorithm used here is based on https://doi.org/10.1122/1.3473925
% i.e. hoda larson journal of rheology 2010
% brownian dynamics ... polymers...
% see pages 1065-1066

try
    load([pfx '_' num2str(irun) '_00sec.mat'])
    load([pfx '_' num2str(irun) '_' num2str(it) 'sec.mat'])
catch
    disp('No file, exiting')
    ncross = nan;
    return
end
rbead0 = rbead;
tact0 = t_rbead;
ifor0 = ifor;

disp('* * * * * * * * * * * Running simulation for 1 step * * * * * * * * * * *')
nstep = 1;
dt = 0.05;
rk1_v48gauss
rbead1 = rbead;
tact1 = t_rbead;
ifor1 = ifor;

ncross = 0;
disp('* * * * * * * * * * * * Searching for crossings * * * * * * * * * * * *')

for ibead = 1:size(rbead0, 2)
    idx_new = find(tact1 == tact0(ibead), 1);
    if isempty(idx_new) || ismember(ibead, ifor0) || ismember(ibead, ifor0+1)
        continue;
    end
    for jbead = ibead+2:size(rbead0, 2)
        jdx_new = find(tact1 == tact0(jbead), 1);
        dij0 = norm(rbead0(:, ibead) - rbead0(:, jbead));
        if(isempty(jdx_new) || ismember(jbead, ifor0) || ismember(jbead, ifor0+1) || dij0 > 2*dbead)
            continue;
        end
        
        % algorithm begins, get all the bead positions
        % and set up some arrays.
        
        d_3 = 0; d_2 = 0; d_1 = 0; d_0 = 0;
        
        x_i = rbead0(1,ibead-1);
        y_i = rbead0(2,ibead-1);
        z_i = rbead0(3,ibead-1);
        
        x_ip1 = rbead0(1,ibead);
        y_ip1 = rbead0(2,ibead);
        z_ip1 = rbead0(3,ibead);
        
        x_j = rbead0(1,jbead-1);
        y_j = rbead0(2,jbead-1);
        z_j = rbead0(3,jbead-1);
        
        x_jp1 = rbead0(1,jbead);
        y_jp1 = rbead0(2,jbead);
        z_jp1 = rbead0(3,jbead);
        
        xn_i = rbead1(1,idx_new-1);
        yn_i = rbead1(2,idx_new-1);
        zn_i = rbead1(3,idx_new-1);
        
        xn_ip1 = rbead1(1,idx_new);
        yn_ip1 = rbead1(2,idx_new);
        zn_ip1 = rbead1(3,idx_new);
        
        xn_j = rbead1(1,jdx_new-1);
        yn_j = rbead1(2,jdx_new-1);
        zn_j = rbead1(3,jdx_new-1);
        
        xn_jp1 = rbead1(1,jdx_new);
        yn_jp1 = rbead1(2,jdx_new);
        zn_jp1 = rbead1(3,jdx_new);
        
        
        x_ip1_i = x_ip1 - x_i;
        y_ip1_i = y_ip1 - y_i;
        z_ip1_i = z_ip1 - z_i;
        
        x_j_i = x_j - x_i;
        y_j_i = y_j - y_i;
        z_j_i = z_j - z_i;
        
        x_jp1_j = x_jp1 - x_j;
        y_jp1_j = y_jp1 - y_j;
        z_jp1_j = z_jp1 - z_j;
        
        d_x_i = xn_i - x_i;
        d_y_i = yn_i - y_i;
        d_z_i = zn_i - z_i;
        
        d_x_j = xn_j - x_j;
        d_y_j = yn_j - y_j;
        d_z_j = zn_j - z_j;
        
        d_x_ip1 = xn_ip1 - x_ip1;
        d_y_ip1 = yn_ip1 - y_ip1;
        d_z_ip1 = zn_ip1 - z_ip1;
        
        d_x_jp1 = xn_jp1 - x_jp1;
        d_y_jp1 = yn_jp1 - y_jp1;
        d_z_jp1 = zn_jp1 - z_jp1;
        
        
        d_x_ip1_i = d_x_ip1 - d_x_i;
        d_y_ip1_i = d_y_ip1 - d_y_i;
        d_z_ip1_i = d_z_ip1 - d_z_i;
        
        d_x_j_i = d_x_j - d_x_i;
        d_y_j_i = d_y_j - d_y_i;
        d_z_j_i = d_z_j - d_z_i;
        
        d_x_jp1_j = d_x_jp1 - d_x_j;
        d_y_jp1_j = d_y_jp1 - d_y_j;
        d_z_jp1_j = d_z_jp1 - d_z_j;
        
        
        % following a_1,b_1,c_1 and the loop
        % are all ways of inputting the massive
        % eq.(12) into our program.
        
        c_1 = [1,1,1,-1,-1,-1];
        
        a_1 = [x_j_i,x_ip1_i,x_jp1_j,x_j_i,x_jp1_j,x_ip1_i];
        a_2 = [y_jp1_j,y_j_i,y_ip1_i,y_ip1_i,y_j_i,y_jp1_j];
        a_3 = [z_ip1_i,z_jp1_j,z_j_i,z_jp1_j,z_ip1_i,z_j_i];
        %a_3 = [z_ip1_i,z_jp1_j,z_jp1_j,z_jp1_j,z_ip1_i,z_j_i];
        % commented a_3 above was the formula in the paper, but I think there is a typo
        
        b_1 = [d_x_j_i,d_x_ip1_i,d_x_jp1_j,d_x_j_i,d_x_jp1_j,d_x_ip1_i];
        b_2 = [d_y_jp1_j,d_y_j_i,d_y_ip1_i,d_y_ip1_i,d_y_j_i,d_y_jp1_j];
        b_3 = [d_z_ip1_i,d_z_jp1_j,d_z_j_i,d_z_jp1_j,d_z_ip1_i,d_z_j_i];
        %b_3 = [d_z_ip1_i,d_z_jp1_j,d_z_jp1_j,d_z_jp1_j,d_z_ip1_i,d_z_j_i];
        % commented b_3 above was the formula in the paper, but I think there is a typo
        
        for k=1:6
            d0 = a_1(k) * a_2(k) * a_3(k);
            
            d1 = a_1(k) * a_2(k) * b_3(k);
            d1 = d1 + a_2(k) * a_3(k) * b_1(k);
            d1 = d1 + a_3(k) * a_1(k) * b_2(k);
            
            d2 = a_1(k) * b_2(k) * b_3(k);
            d2 = d2 + a_2(k) * b_3(k) * b_1(k);
            d2 = d2 + a_3(k) * b_1(k) * b_2(k);
            
            d3 = b_1(k) * b_2(k) * b_3(k);
            
            d_3 = d_3 + c_1(k) * d3;
            d_2 = d_2 + c_1(k) * d2;
            d_1 = d_1 + c_1(k) * d1;
            d_0 = d_0 + c_1(k) * d0;
        end
        
        bb = roots([d_3,d_2,d_1,d_0]); % this is equation (13) of the paper.
        alp = 0;
        
        for cnt2 = 1:length(bb)
            rts = bb(cnt2);
            if isreal(rts)
                if rts > 0 & rts < 1
                    alp = rts;
                end
            end
        end
        
        % following 8 lines or so are eq(12) of the paper
        x_ip1 = x_ip1 + alp*(xn_ip1 - x_ip1);
        y_ip1 = y_ip1 + alp*(yn_ip1 - y_ip1);
        x_jp1 = x_jp1 + alp*(xn_jp1 - x_jp1);
        y_jp1 = y_jp1 + alp*(yn_jp1 - y_jp1);
        
        x_i = x_i + alp*(xn_i - x_i);
        y_i = y_i + alp*(yn_i - y_i);
        x_j = x_j + alp*(xn_j - x_j);
        y_j = y_j + alp*(yn_j - y_j);
        
        x_ip1_i = x_ip1 - x_i;
        y_ip1_i = y_ip1 - y_i;
        x_j_i = x_j - x_i;
        y_j_i = y_j - y_i;
        x_jp1_j = x_jp1 - x_j;
        y_jp1_j = y_jp1 - y_j;
        
        d_r = x_ip1_i * y_jp1_j - y_ip1_i * x_jp1_j;
        n_r_1 = x_j_i * y_jp1_j - y_j_i * x_jp1_j;
        n_r_2 = x_j_i * y_ip1_i - y_j_i * x_ip1_i;
        
        % eqs. (10) & (11) of the paper.
        p_i = n_r_1 / d_r;
        p_j = n_r_2 / d_r;
        
        if (alp > 0 & alp < 1) & (p_i > 0 & p_i < 1) & ...
                (p_j > 0 & p_j < 1)
            ncross = ncross + 1;
            
            figure
            hold on
            p = [rbead0(:, ibead-1) rbead0(:, ibead)];
            q = [rbead0(:, jbead-1) rbead0(:, jbead)];
            plot3(p(1, :), p(2, :), p(3, :), 'ro-', 'linewidth', 2)
            plot3(q(1, :), q(2, :), q(3, :), 'bo-', 'linewidth', 2)
            p = [rbead1(:, idx_new-1) rbead1(:, idx_new)];
            q = [rbead1(:, jdx_new-1) rbead1(:, jdx_new)];
            plot3(p(1, :), p(2, :), p(3, :), 'ro--', 'linewidth', 2)
            plot3(q(1, :), q(2, :), q(3, :), 'bo--', 'linewidth', 2)
            grid on
            view([5 5 5])
            axis equal
            close;
            
%             return;
            %}
        end
    end
end
