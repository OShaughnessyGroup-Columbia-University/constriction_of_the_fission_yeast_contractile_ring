function [vm] = vnode(diry)
	% find all '_10min.mat' files in directory diry
	cc = dir([diry '\' '*10min.mat']);
	vm = [];
	
	for kk=1:length(cc)
	
		% extract velocities and myosin positions;
		load([cc(kk).folder '\' cc(kk).name],...
			'vmyo','rmyo',...
			'rbead','ifor','ipt',...
			'fmmat','r_ring');
			
		rbead(3,:) = rbead(3,:) - mean(rbead(3,:));
		rmyo(3,:) = rmyo(3,:) - mean(rmyo(3,:));
		
		[theta,rho,z] = cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
		theta = theta + pi/2;
		dirn = [cos(theta);sin(theta);0*(theta)];
		vmtemp = (sum(vmyo .* dirn,1));
		
		% find those myosins that just have one formin
		% attached to them.
		nfor_per_myo = sum(fmmat,1);
		uni_for_per_myo = find(nfor_per_myo == 1);
		
		% find these formins
		[for_of_myo,junk] = find(fmmat(:,uni_for_per_myo));
		
		% loop over these filaments and see if they are suitable
		good_flt = true(1,length(for_of_myo));
		cnt = 1;
		
		for ll=1:length(for_of_myo)
			curr_flt = rbead(:,ifor(for_of_myo(ll)):...
									ipt(for_of_myo(ll)));
			if(any(abs(curr_flt(3,:)) > 100/1000))
				good_flt(cnt) = false;
			end
			if(any(sqrt(curr_flt(1,:).^2 + curr_flt(2,:).^2) < r_ring - 200/1000))
				good_flt(cnt) = false;
			end
			if(size(curr_flt,2) < 3)
				good_flt(cnt) = false;
			else
				x1 = curr_flt(1,2);
				y1 = curr_flt(2,2);
				x2 = curr_flt(1,3);
				y2 = curr_flt(2,3);
				v1 = [x2-x1;y2-y1;];
				v2 = [x2+x1;y2+y1;];
				if(sum(v1.*v2) > 0)
					good_flt(cnt) = false;
				end
			end
			cnt = cnt + 1;
		end
		
		temp1 = vmtemp(uni_for_per_myo);
		vm = [vm temp1(good_flt)];
		
	end
end