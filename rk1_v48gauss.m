    istore = 1:nstep;
    state = zeros(size(istore));
%     state = cell(length(istore),11);
if ~exist('bmmat','var')
    bpmat_trial = false(size(rbead,2),size(rmyo,2));
    bmmat = bm_nbr (rbead,rmyo,ifor,ipt,rcapmyo_long,rcapmyo_short,fmmat);
end
if ~exist('bpmat','var')
    bpmat = bp_nbr (rbead,rmyp,ifor,ipt,rcapmyp);
end

    for itimes = 1:nstep
% %         % update state parameters other than positions
% %         poly_cycle = round(dbead/vpol/dt);  % number of time steps for an actin filament to add one bead
% %         if mod(itimes,poly_cycle) == 0
% %             poly_time = true;
% %         else
% %             poly_time = false;
% %         end
        [rbead,rmyo,rmyp,xmat,fmmat,ipt,ifor,bancf,bancm,dbead_first,bpmat_trial,t_rbead,t_rmyo,t_rmyp] = update_v23(rbead,rmyo,rmyp,bpmat,xmat,rhoain,rxbind,rex_act,rexc_n,fmmat,bancf,bancm,ifor,ipt,kofffor,...
            koffmyo,koffmyp,dt,d_for,d_myo,dbead,dbead_first,rsev,koffx,rcapmyp,l_break_sq,r_ring,wr,...
            binding_rate_myo2,binding_rate_myp2,binding_rate_for,vpol,geometry,t,t_rbead,t_rmyo,t_rmyp,...
			alphafor,epsilonfor,qmax);
        bmmat = bm_nbr (rbead,rmyo,ifor,ipt,rcapmyo_long,rcapmyo_short,fmmat);
        bpmat = bp_nbr (rbead,rmyp,ifor,ipt,rcapmyp);
        rdiff = get_rdiff(rbead,ifor,ipt); 
        nsat_p = nsat;
        [force,force_vec,gbeadfv,gppf,bmmatpt,bpmatpt]=get_force_v15gauss(rbead,rmyo,rmyp,rdiff,ifor,ipt,bancm,bmmat,bpmat,nhead,fhead,fheadmyp,fone,nsat,rcapmyo_u,rcapmyo_long,rcapmyo_short,rcapmyp,kcap,kcap_p,r_flat_myo,r_flat_myp,dbead,vmyo0,kod2,r_ring,kwall,kexc_oo,kexc_pp,kexc_po,rexc,rexc_n,rexc_po,kex_act,rex_act,xmat,rx0,kx,d_cdc15,d_myo,geometry);
        amat = get_a_v1(rbead,rmyo,rmyp,bmmat,bpmat,gppf,rdiff,dt,ga_myo2,ga_myp2,ifor,ipt,bmmatpt,bpmatpt);
        if geometry == 'cyl'
            r_ring = r_ring - vs_abs * dt;
            [pcpq,cvec,cdot] = get_pcpq_xy(rbead,rmyo,rmyp,ifor,fmmat,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo);
        elseif geometry == 'sph'
            zbm = mean(rbead(3,:));
            r_ring = sqrt(r*r - zbm*zbm);
            [pcpq,cvec,cdot] = get_pcpq_xy(rbead,rmyo,rmyp,ifor,fmmat,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo,d_stem);
        end        
%         w = get_w(ga,ifor,ipt,rbead,rmyo,bancf,bancm,bmmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt);
        m = get_m_v1 (dt,ga_myo2,ga_myp2,ifor,ipt,rbead,rmyo,rmyp,bancf,bancm,bmmat,bpmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt,bpmatpt);
        m = m - amat;
        nbead = size(rbead,2);
        nanc = sum(bancf)+sum(bancm);
        anc_flag = false(1,size(pcpq,1));
        anc_flag(1:nanc) = true;
        [vbead,vmyo,vmyp,fc,fanc] = velocity(pcpq,m,force,nbead,rmyo,cvec,cdot,tau,anc_flag);
        % evolve
        rbead = rbead + dt * vbead;
        rmyo = rmyo + dt * vmyo;
        rmyp = rmyp + dt * vmyp;
		
		if debug_cnt >= 0
			save(strcat('dump/',prefix,num2str(itask),'_',num2str(debug_cnt),'.mat'), ...
			'rbead','rmyo','rmyp','ifor','fmmat',...
			'r','r_ring','fc','fanc',...
			't_rbead','t_rmyo','t_rmyp','dt');
			debug_cnt = debug_cnt + 1;
		end
    end
