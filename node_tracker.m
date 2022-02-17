function node_tracker(irun, pfx, dr, tlist)
%irun = 3;
%tlist = 5:5:120;
%pfx = 'kxscan4';
%load([pfx '_' num2str(irun) '_00sec.mat']);
status = mkdir('images');
load([dr '/' pfx '_' num2str(irun) '_' num2str(tlist(1)) 'sec.mat']);
t_init = t_rmyo;
ipt = [ifor(2:end)-1 size(rbead,2)];
dt = 1;
% t0 = t_rmyo(randi(length(rmyo))); % choose tag of random node.

for i0 = 1:length(rmyo)
    t0 = t_init(i0); % choose tag of random node.
    tbreak = max(tlist);
    for it = tlist
	imgName = ['./images/' pfx '_' num2str(irun) '_node' num2str(i0) '_' num2str(it) 'sec.png'];
	if isfile(imgName)
	    continue;
	end
        load([dr '/' pfx '_' num2str(irun) '_' num2str(it) 'sec.mat']);
	ipt = [ifor(2:end)-1 size(rbead,2)];
        inode = find(t_rmyo == t0);
        if(isempty(inode))
            tbreak = it;
            break
        end
        fig = figure;
        xlim([-2, 2])
        ylim([-2, 2])
        axis square
        hold on
        plot(r_ring*cos(0:0.01:2*pi), r_ring*sin(0:0.01:2*pi), 'k')
        plot(rmyo(1, inode), rmyo(2, inode), 'ro', 'Markersize', 18);
        % plot3(rmyo(1, inode), rmyo(2, inode), rmyo(3, inode), 'ro', 'Markersize', 18);
        flist = find(fmmat(:, inode));
        for idx_fil = flist'
            plot(rbead(1, ifor(idx_fil):ipt(idx_fil)),...
                  rbead(2, ifor(idx_fil):ipt(idx_fil)), '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
            %plot3(rbead(1, ifor(idx_fil):ipt(idx_fil)),...
            %    rbead(2, ifor(idx_fil):ipt(idx_fil)), ...
            %    rbead(3, ifor(idx_fil):ipt(idx_fil)), 'o-');
        end
	hold on;
	%text(0.0,0.0,[num2str(dt*it) 's'])
        print(fig,imgName,'-dpng');
        close
    end
    vidName = ['./images/' pfx '_' num2str(irun) '_node' num2str(i0) '.avi'];
    if isfile(vidName)
	continue;
    end
    vid = VideoWriter(vidName);
    vid.FrameRate = 10;
    % vid.FrameRate = floor(1/dt);
    open(vid);
	for it = tlist(tlist < tbreak);
	   img = imread(['./images/' pfx '_' num2str(irun) '_node' num2str(i0) '_' num2str(it) 'sec.png']);
	   writeVideo(vid,img)
    end
    close(vid)
end
end
