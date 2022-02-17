% clear
% close
% %% load parameters from saved initial file
% % import from steady state configuration of Matt's simulation
%     name = 'maxfb_4s'; % tail can be randomized if needed
%     load(name);

%load maxfb_0ds
%% RK1 (Euler method) to evolve r and v

    istore = 1:nstep;
    state = zeros(size(istore));
%     state = cell(length(istore),11);
    for itimes = 1:nstep
        itimes
        rdiff = get_rdiff(rbead,ifor,ipt);
        bmmat = bm_nbr (rbead, rmyo, ipt, rcap);
        for ifv = 1:1e6
            [force, gbeadfv, gppf, bmmatpt] = get_force(maxfb,rbead,rmyo,rdiff,ifor,ipt,bmmat,nhead,fhead,fone,nsat,kcap,vmyo0,kod2,lr,kwall,kexc,rexc,xmat);
            amat = get_a(rbead, rmyo, bmmat, gppf, rdiff, ga, ifor, ipt, bmmatpt);
            [pcpq, cvec] = get_pcpq(rbead,rmyo,ifor,bancf,bancm,xmat,r_ring_sq,dbead_sq);
            w = get_w (ga,ifor,ipt,rbead,rmyo,bancf,bancm,bmmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt);
            [vbead, vmyo,fc] = velocity(pcpq,w,force,nbead,amat,alpha, cvec);
            bmmat_new = modfv(bmmat, vbead, vmyo, rdiff, vmyo0);        % when the relative velocity between a myosin and an actin filament exceeds vmyo0, there is no pulling force between them
            if ~any(any(bmmat_new ~= bmmat))    % if they are equal
                break
            else
                bmmat = bmmat_new;
            end
        end
        % evolve
        rbead = rbead + dt * vbead;
        rmyo = rmyo + dt * vmyo;
        % update state parameters other than positions
        [vbead,rbead,vmyo,rmyo,xmat,ipt,ifor,nbead,bancf,bancm] = update (rbead,rmyo,vbead,vmyo,xmat,bancf,bancm,ifor,ipt,kofffor,koffmyo,dt,dbead,rsev,koffx,rxbind,fc,fpeel,l_break_sq);    
    end