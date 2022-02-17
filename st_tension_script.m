rdiff = get_rdiff(rbead,ifor,ipt);
        bmmat = bm_nbr (rbead, rmyo, ifor, ipt, rcapmyo_long, rcapmyo_short, fmmat);
        bpmat = bp_nbr (rbead, rmyp, ifor, ipt, rcapmyp);
        nsat_p = nsat;
        [force, gbeadfv, gppf, bmmatpt, bpmatpt] = get_force(rbead,rmyo,rmyp,rdiff,ifor,ipt,bancm,bmmat,bpmat,nhead,fhead,fheadmyp,fone,nsat,rcapmyo_u, rcapmyo_long, rcapmyo_short,rcapmyp,kcap,kcap_p,dbead,vmyo0,kod2,r_ring,kwall,kexc,rexc,xmat,rx0,kx,geometry);
        amat = get_a(rbead, rmyo, rmyp, bmmat, bpmat, gppf, rdiff, ga, ifor, ipt, bmmatpt, bpmatpt);
                    [pcpq, cvec, cdot] = get_pcpq_xy(rbead,rmyo,rmyp,ifor,fmmat,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo);
m = get_m (ga,ifor,ipt,rbead,rmyo,rmyp,bancf,bancm,bmmat,bpmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt,bpmatpt);
        m = m - amat;
        tau = tau_o_dt * 0.01;
        nbead = size(rbead,2);
        nanc = sum(bancf)+sum(bancm);
        anc_flag = false(1,size(pcpq,1));
        anc_flag(1:nanc) = true;
[vbead, vmyo, vmyp, fc_sept] = ...
	st_velocity(pcpq,m,force,nbead,size(rmyo,2),cvec,cdot,tau,anc_flag);
[vbead, vmyo, vmyp, fc] = velocity(pcpq,m,force,nbead,rmyo,cvec,cdot,tau);
fanc = vec2mat(fc_sept,3)';
fanc_j = vec2mat(fc,3)';
% calculate roughness, approx
					
            % calculate tension, approx method
				ni = -[rbead(:,ifor),rmyo];
				ni = ni./repmat(sqrt(sum(ni.*ni,1)),[3 1]);
				ni(3,:) = 0;
				
				% some files have a different fanc structure
				% if so, correct
				if(size(fanc,2) > size(ni,2))
					fanc = fanc(:,[bancf,bancm]);
                    fanc_j = fanc_j(:,[bancf,bancm]);
					% fanc is 3 by N array equivalent to fc_sept
					% in velocity.m
				end
				
				% calculate tension and append
				temp1 = 1/(2*pi) * sum(dot(-fanc,ni));
				tens_v = temp1
