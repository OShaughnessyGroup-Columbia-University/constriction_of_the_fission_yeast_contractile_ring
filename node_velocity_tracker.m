function node_velocity_tracker(irun, pfx,dr,tlist)
	status = mkdir('plots');

	load([dr '/' pfx '_' num2str(irun) '_' num2str(tlist(1)) '.mat']);
	t_rmyo0 = t_rmyo;
	rmyo0 = rmyo;

	for i0 = 1:length(rmyo)
	    v = nan(1,length(tlist));
	    rmyoOld = nan(3,1);

	    t0 = t_rmyo(i0); % choose tag of random node.
	    cnt = 1;

	    for it = tlist
		load([dr '/' pfx '_' num2str(irun) '_' num2str(it) '.mat']);
		inode = find(t_rmyo == t0);
		if(isempty(inode))
		    break
		end
		v(cnt) = norm(rmyo(:,inode) - rmyoOld)/dt * 1000;

		if cnt > 1
		    thetaHat = [rmyo(2,inode),-rmyo(1,inode),0]';
		    if sum(thetaHat.*(rmyo(:,inode)-rmyoOld)) < 0
			v(cnt) = -v(cnt);
		    end
		end

		cnt = cnt + 1;
		rmyoOld = rmyo(:,inode);
            end

	    fig = figure;
	    okIndx = ~isnan(v);
	    % need at least 10 points in the v vs t curve
	    if length(okIndx) < 10
		continue;
	    end
	    plot(dt*tlist(okIndx),v(okIndx),'b.-');
	    xlim([dt*min(tlist(okIndx)) dt*max(tlist(okIndx))]);
	    xlabel('Time (s)');
	    ylabel('Myo2 node velocity (nm/s)');
	    savefig(['./plots/' pfx '_' num2str(irun) '_node' num2str(i0) '_v.fig']);
	    close
	end
end
