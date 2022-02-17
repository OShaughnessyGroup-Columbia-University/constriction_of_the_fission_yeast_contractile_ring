function [] = st_make_2d_movie(itask,flag)

	% scan all files and obtain the amplitude
	% of the sinusoidal distortion that is slowly
	% corrected by the mechanosensitivity of the Bgs

	% some file information and setting up arrays
	prefix0 = 'noga';
	N = 1200;
	
	% image names tension and time
	imageNames = cell(1,N);
% 	Tarray = zeros(1,N);
 	tarray = zeros(1,N);
% 	Rarray = zeros(1,N);
% 	wLarray = zeros(1,N);
	
	namecnt = 1;
	
	%N = 549;
% 	amp = zeros(1,N);
	% amplitude of distortion
% 	rmean = zeros(1,N);
	% mean radius of septumedge

	workingDir = '.';
	imageDir = [workingDir '/' 'images'];
	movieDir = [workingDir '/' 'movies'];
	plotDir = [workingDir '/' 'plots'];

	if(exist(imageDir)~=7)
		mkdir(imageDir)
	end

	if(exist(movieDir)~=7)
		mkdir(movieDir)
	end

	if(exist(plotDir)~=7)
		mkdir(plotDir)
	end
	
	% load latest result
	mat_files = dir([workingDir '/' sprintf('%s_%d_*%s.mat',prefix0,itask,'sec')]);
	fname_sel = cell(1,length(mat_files));
	
	for k=1:length(mat_files)
		fname_sel{k} = mat_files(k).name;
	end
% 	mat_files
	fname_sel = natsortfiles(fname_sel);

	% loop over files
	for k2 = 1:length(mat_files)
%     for k2 = 1:90
		fnamemat = fname_sel{k2};
		imageNames{namecnt} = sprintf('%sF.png',fnamemat(1:(length(fnamemat)-4)));
		if(flag)
			imageNames{namecnt} = sprintf('%s.png',fnamemat(1:(length(fnamemat)-4)));
		end
	
		if (exist(fullfile(imageDir,imageNames{namecnt}), 'file') == 2)
			namecnt = namecnt + 1;
			continue;
		else
			load(fnamemat);
			disp(fnamemat);
		end
		
		% below line necessary as matlab confuses
		% internal fn xmat with variable xmat
		load(fnamemat,'xmat');
		
		% combine anchored flag of formin and myosin
		bancmrd = zeros(1,size(bancf,2)+size(bancm,2));
		bancmrd(1,length(bancf)+1:end) = bancm;
		bancmrd = logical(bancmrd);

		% extract positions of formins
		rfor = rbead(:,bancf);
		
		% extract crosslinker locations
		[xrow,xcol,xval] = find(tril(xmat));

		% plot anc myo, form and the reference ring
		% w.r.t. arrows 1nm = 1pN
		% remember all forces in the plot are forces from the membrane on the ring anchors
		% so on average they have to be outward.

		%pause(0.1);
		fig=figure('visible','off');
		%clf;
		nfor = length(ifor);
		lactmn = 0.1 * (size(rbead,2)-nfor)/nfor;
		
		if(flag)
			for k = 1:nfor
				plot3(rbead(1,ifor(k):ipt(k)),rbead(2,ifor(k):ipt(k)),rbead(3,ifor(k):ipt(k)),'b-','LineWidth',1);
				hold on;
			end
			plot3(rbead(1,ifor(k):end),rbead(2,ifor(k):end),rbead(3,ifor(k):end),'b-','LineWidth',1);
			plot3(rfor(1,:),rfor(2,:),rfor(3,:),'g.');
		end
		
		nmyoa = sum(bancm);
		nmyou = size(rmyp,2)/2;
		
		if(flag)
			plot3(rmyo(1,bancm),rmyo(2,bancm),rmyo(3,bancm),'ro');
			hold on;
% 			for k3=1:2:size(rmyp,2)
			for k3=1:size(rmyp,2)
				plot3(rmyp(1,k3),rmyp(2,k3),rmyp(3,k3),'r*-');
% 				plot3(rmyp(1,k3:k3+1),rmyp(2,k3:k3+1),rmyp(3,k3:k3+1),'r*-');
			end
			hold on;
		end
		
			
		view([0 90]);
		
		% plot crosslinkers
		if(flag)
			nx = length(xrow);
			for k = 1:nx
				plot3(rbead(1,[xrow(k) xcol(k)]),rbead(2,[xrow(k) xcol(k)]),rbead(3,[xrow(k) xcol(k)]),'k-','LineWidth',1);
				hold on;
			end
		end
		
		%Virial tension
		
		%{
		% wrong tension calculation
		indx = nearestNeighbor(T0,p_sept');
		ni = -vertexNormal(T0,indx)';
		ni(3,:) = 0;
		T = 1/(2*pi) * sum(dot(-fanc,ni));
		%}
		
		axis equal;
		xlim([-2.5 2.5]);
		ylim([-2.5 2.5]);
		xlabel('x (um)');
		ylabel('y (um)');
		
		%Tarray(namecnt) = T;
% 		tarray(namecnt) = t;
		%{
		Rarray(namecnt) = R;
		wLarray(namecnt) = wL;
		%}
		
		%leg_keys = {'t','T','L','nfor','nmyoa','nmyou','lactmn','nx','R','wL','xc','yc'};
		leg_keys = {'nfor','nmyoa','nmyou','lactmn','nx'};
		
		% pad strings in leg_keys properly
		
% 		a = strvcat(leg_keys{:});
% 		leg_keys = mat2cell(a,ones(size(a,1),1),size(a,2));
		%{
		t1 = char(leg_keys);
		for k=1:length(leg_keys)
			leg_keys{k} = t1(k,:);
		end
		%}
		
		%leg_values = [k2,round(T),nfor,nmyoa,nmyou,lactmn,nx,R,wL*1000,xc*1000,yc*1000];
		%leg_units = {'s','pN','um','','','','um','','um','nm','nm','nm'};
		%leg_format = {'%d','%d','%.1f','%d','%d','%d','%.1f','%d','%.2f','%.2f','%.2f','%.2f'};
		
		if(flag)
% 			leg_values = [t,nfor,nmyoa,nmyou,lactmn,nx];
% 			leg_units = {'s','','','','um',''};
% 			leg_format = {'%.2f','%d','%d','%d','%.1f','%d'};
			leg_values = [nfor,nmyoa,nmyou,lactmn,nx];
			leg_units = {'','','','um',''};
			leg_format = {'%d','%d','%d','%.1f','%d'};
			
			leg_str = cell(1,length(leg_units));
			
			for k=1:length(leg_str)
				leg_str{k} = [leg_keys{k} ' ' num2str(leg_values(k),leg_format{k}) ' ' leg_units{k}];
			end
			
			annotation('textbox',...
				[0.78 0.45 0.22 0.4],...
				'String',leg_str,...
				'FontSize',8,...
				'FontName','Courier',...
				'LineStyle','-',...
				'EdgeColor',[0 0 0],...
				'LineWidth',2,...
				'BackgroundColor',[0.9  0.9 0.9],...
				'Color',[0.84 0.16 0]);
		end
			
		print(fig,fullfile(imageDir,imageNames{namecnt}),'-dpng');
		
		clf;
		namecnt = namecnt + 1;
		
	end
	hold off;
	close all;
	
	namecnt = namecnt - 1;
	
	%Tarray = Tarray(1:namecnt);
	tarray = tarray(1:namecnt);
	
	%{
	Rarray = Rarray(1:namecnt);
	wLarray = wLarray(1:namecnt);
	%}

	% make the movie

	% code mostly copied from 
	% https://www.mathworks.com/help/matlab/examples/convert-between-image-sequences-and-video.html

	%Create New Video with the Image Sequence
	%Construct a VideoWriter object, which creates a Motion-JPEG AVI file by default.

	outputVideo = [];
	
	if(flag)
		outputVideo = VideoWriter(fullfile(movieDir,sprintf('%s_%d.avi',prefix0,itask)));
	else
		outputVideo = VideoWriter(fullfile(movieDir,sprintf('%s_%d%s.avi',prefix0,itask,'F')));
	end
	
	outputVideo.FrameRate = 5;
	open(outputVideo)

	%Loop through the image sequence, load each image, and then write it to the video.

	for ii = 1:namecnt
	   img = imread(fullfile(imageDir,imageNames{ii}));
	   writeVideo(outputVideo,img)
	end

	%Finalize the video file.
	close(outputVideo)
	
	close all;
	%{
	% plot tension versus time
	fig=figure('visible','off');
	plot(tarray,Tarray,'bo-');
	xlabel('t (s)');
	ylabel('Tension (pN)');
	print(fig,fullfile(plotDir,sprintf('%s_%d_tens.png',prefix,itask)),'-dpng');
	clf;
	
	% plot ring radius versus time
	plot(tarray,Rarray,'bo-');
	xlabel('t (s)');
	ylabel('R (um)');
	print(fig,fullfile(plotDir,sprintf('%s_%d_R.png',prefix,itask)),'-dpng');
	clf;
	% plot 'roughness' versus time
	plot(tarray,wLarray*1000,'bo-');
	xlabel('t (s)');
	ylabel('Roughness w_L (nm)');
	print(fig,fullfile(plotDir,sprintf('%s_%d_wL.png',prefix,itask)),'-dpng');
	clf;
	save(fullfile(plotDir,sprintf('%s_%d_data.mat',prefix,itask)));
	%}
	
	close all;
	
end
