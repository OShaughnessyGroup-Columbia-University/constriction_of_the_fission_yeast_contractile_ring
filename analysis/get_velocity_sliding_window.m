function v = get_velocity_sliding_window(prefix,twindow,tstart,tend,tstep,runIndx);

	% load files with prefix and runIndx indices given by the eponymous
	% array. within these files, construct windows of t to t+twindow
	% and measure node motions. slide t in steps of tstep from
	% tstart to tend.
	
	v = [];
	for k=runIndx
		for t1 = tstart:tstep:tend
			f1 = ['coopf_' num2str(k) '_' num2str(t1) 'sec.mat'];
			f2 = ['coopf_' num2str(k) '_' num2str(t1+twindow) 'sec.mat'];
			if isfile(f1) & isfile(f2)
				v2 = get_velocity_bw_files(f1,f2);
				v = [v v2];
			end
		end
	end
end