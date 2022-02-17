function v = get_velocity_bw_files(f1,f2)
	% given two filenames f1, f2
	% find cdc12 nodes that persisted between these files
	% and get the magnitude of their velocities.
	
	% load the data
	
	[t_rbead1,rbead1,vbead,t1] = load_data(f1);
	
	if strcmp(f1,f2)
		v = sqrt(sum(vbead.*vbead,1));
		return;
	end
	
	[t_rbead2,rbead2,vbead,t2] = load_data(f2);
	
	[B,C] = ismember(t_rbead1,t_rbead2);
	dsp = nan(1,length(t_rbead1));
	
	for k=1:length(t_rbead1)
		if C(k) > 0
			dr = (rbead1(:,k)-rbead2(:,C(k)));
			dsp(k) = norm(dr);
		end
	end
	
	dt = t2 - t1;
	
	v = dsp(~isnan(dsp))./dt;
	
end

function [t_rbead1,rbead1,vbead,tstamp] = load_data(f1);

	load(f1,'t_rbead','ifor','rbead','vbead','fmmat');
	tstamp = max(t_rbead);
	t_rbead1 = t_rbead(ifor);
	rbead1 = rbead(:,ifor);
	vbead = vbead(:,ifor);
	
	% 1 node may have many formins. in that case, measure
	% the velocity of only one of those formins.
	forFlag = false(1,length(ifor));
	for k=1:size(fmmat,2)
		B = find(fmmat(:,k));
		if length(B) > 0
			B = B(randperm(length(B)));
			forFlag(B(1)) = true;
		end
	end
	
	t_rbead1 = t_rbead1(forFlag);
	rbead1 = rbead1(:,forFlag);
	vbead = vbead(:,forFlag);
end