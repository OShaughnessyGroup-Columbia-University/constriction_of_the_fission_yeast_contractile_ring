function [] = st_calculate_tension(dirName,itask)

	currDir = pwd;
	cd(dirName);
	% load file

	dir_list = {''};

	root_dir = '.'

	prefix = 'test';
	suffix = 'us';
	f_per_unit_len = [];
	curv = [];
	ang_f_curv = [];
	tens_v = [];
	t_v = [];
	wL_v = [];

	workingDir = '.';
	imageDir = [workingDir '/' 'images'];
	movieDir = [workingDir '/' 'movies'];
	plotDir = [workingDir '/' 'plots'];
	
	prefix0 = 'test';
	
	if(exist(imageDir)~=7)
		mkdir(imageDir);
	end

	if(exist(movieDir)~=7)
		mkdir(movieDir);
	end

	if(exist(plotDir)~=7)
		mkdir(plotDir);
	end

	for j = 1:length(dir_list)

		result_dir = [root_dir dir_list{j} '/'];

		for runnum = itask

			% load latest result
			mat_files = dir([result_dir sprintf('*_%d_*%s.mat',runnum,suffix)]);
			fname_sel = cell(1,length(mat_files));
			
			for k=1:length(mat_files)
				fname_sel{k} = mat_files(k).name;
			end
			
			fname_sel = natsortfiles(fname_sel);
			
			namecnt = 1;
			
			for k=1:length(mat_files)
				load([result_dir fname_sel{k}]);
					
				T0=triangulation(t0,p0);
				% calculate roughness, approx
					
				% calculate tension, approx method
				ni = -[rbead(:,ifor),rmyo];
				ni = ni./repmat(sqrt(sum(ni.*ni,1)),[3 1]);
				%ni(3,:) = 0;
				
				% some files have a different fanc structure
				% if so, correct
				if(size(fanc,2) > size(ni,2))
					fanc = fanc(:,[bancf,bancm]);
					% fanc is 3 by N array equivalent to fc_sept
					% in velocity.m
				end
				
				% calculate tension and append
				temp1 = 1/(2*pi) * sum(dot(-fanc,ni));
				tens_v = [tens_v temp1];
				t_v = [t_v t];
				disp(sprintf('%d',k))
			end
			close all;
			save(fullfile(plotDir,sprintf('%s_%d_data.mat',prefix,runnum)));
			fig = figure(1);
			plot(t_v(tens_v >0 & tens_v<1000),tens_v(tens_v >0 & tens_v<1000),'bo-','LineWidth',1);
			hold on;
			xlabel('t (s)');
			ylabel('T (pN, approx, cut outl.)');
			print(fig,fullfile(plotDir,sprintf('%s_%d_tension.png',prefix,runnum)),'-dpng');
			clf;
			
			namecnt = namecnt - 1;
			hold off;
			close all;
			
		end
			
	end
	cd(currDir);
end