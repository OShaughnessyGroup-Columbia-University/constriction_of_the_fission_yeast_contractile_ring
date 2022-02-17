%% RK1 (Euler method) to evolve r and v

    istore = 1:nstep;
    state = zeros(size(istore));
%     state = cell(length(istore),11);
if ~exist('bmmat','var')
    bmmat_trial = false(size(rbead,2),size(rmyo,2));
    bmmat = bm_nbr (rbead, rmyo, ipt, rcap, bancm, select4, atp_cycle, dt, bmmat_trial);
end
    for itimes = 1:nstep
        itimes
% % %         % update state parameters other than positions
% % %         poly_cycle = round(dbead/vpol/dt);  % number of time steps for an actin filament to add one bead
% % %         if mod(itimes,poly_cycle) == 0
% % %             poly_time = true;
% % %         else
% % %             poly_time = false;
% % %         end
        %[vbead,rbead,vmyo,rmyo,xmat,ipt,ifor,nbead,bancf,bancm] = update_ver1(rbead,rmyo,vbead,vmyo,xmat,bancf,bancm,ifor,ipt,kofffor,koffmyo,dt,dbead,rsev,koffx,rxbind,fc,fpeel,l_break_sq,r,wr,binding_rate_myo,poly_time);
        [rbead,rmyo,xmat,ipt,ifor,bancf,bancm,dbead_first] = update_xy(rbead,rmyo,vbead,vmyo,xmat,bancf,bancm,ifor,ipt,kofffor,...
    koffmyo,dt,d_for,d_myo,dbead,dbead_first,rsev,koffx,rxbind,rcap,l_break_sq,r_ring,wr,binding_rate_myo2,binding_rate_myp2,binding_rate_for,vpol,binding_rate_x);
        rdiff = get_rdiff(rbead,ifor,ipt);
        bmmat = bm_nbr (rbead, rmyo, ipt, rcap, bancm, select4, atp_cycle, dt, bmmat_trial);
        [force, gbeadfv, gppf, bmmatpt] = get_force(maxfb,rbead,rmyo,rdiff,ifor,ipt,bmmat,nhead,fhead,fone,nsat,kcap,dbead,vmyo0,kod2,r,...
            kwall,kexc,rexc,xmat,rx0,kx);        
        amat = get_a(rbead, rmyo, bmmat, gppf, rdiff, ga, ifor, ipt, bmmatpt);
        r_ring = r_ring - vs_abs * dt;
        [pcpq, cvec, cdot] = get_pcpq_xy(rbead,rmyo,ifor,bancf,bancm,xmat,r_ring,dbead_sq,vs_abs,dbead_first,vpol);
        w = get_w(ga,ifor,ipt,rbead,rmyo,bancf,bancm,bmmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt);
        [vbead, vmyo,fc] = velocity(pcpq,w,force,size(rbead,2),amat,alpha,cvec,cdot);
        % evolve
        rbead = rbead + dt * vbead;
        rmyo = rmyo + dt * vmyo;
    end