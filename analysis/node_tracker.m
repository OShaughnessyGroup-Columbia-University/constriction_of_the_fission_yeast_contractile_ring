irun = 40;
tlist = 120:2:600;
load(['noagm2_' num2str(irun) '_120sec.mat']);
% t0 = t_rmyo(randi(length(rmyo))); % choose tag of random node.

for i0 = [12 25 29 37 43 76 78 86 91 92 96 102 103]
    t0 = t_rmyo(i0); % choose tag of random node.
    tbreak = max(tlist);
    for it = tlist
        load(['noagm2_' num2str(irun) '_' num2str(it) 'sec.mat']);
        inode = find(t_rmyo == t0);
        flist = find(fmmat(:, inode))';
        if(isempty(inode) || isempty(flist))
            tbreak = it;
            break
        end
        fig = figure;
        xlim([-2, 2])
        ylim([-2, 2])
        axis square
        hold on
        plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi), 'k')
        plot3(rmyo(1, inode), rmyo(2, inode), rmyo(3, inode), 'ro', 'Markersize', 18);
        flist = find(fmmat(:, inode))';
        for idx_fil = flist
            plot3(rbead(1, ifor(idx_fil):ipt(idx_fil)),...
                  rbead(2, ifor(idx_fil):ipt(idx_fil)), ...
                  rbead(3, ifor(idx_fil):ipt(idx_fil)), 'o-');
        end
        print(fig,['./images/noagm2_' num2str(irun) '_node' num2str(i0) '_' num2str(it) 'sec.png'],'-dpng');
        close
    end
    vid = VideoWriter(['./images/noagm2_' num2str(irun) '_node' num2str(i0) '.avi']);
    vid.FrameRate = 5;
    open(vid);
	for it = tlist(tlist < tbreak)
	   img = imread(['./images/noagm2_' num2str(irun) '_node' num2str(i0) '_' num2str(it) 'sec.png']);
	   writeVideo(vid,img)
    end
    close(vid)
end