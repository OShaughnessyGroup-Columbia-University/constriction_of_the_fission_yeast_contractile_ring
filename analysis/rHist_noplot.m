% attempt to define the bundle.

% first, remove formins.
rbeadCutFormin = rbead;
rbeadCutFormin(:,ifor) = [];

nxfil = length(rbead)*dbead/2/pi/r_ring;

% then, create histogram of distances from the median actin bead location in the cross-section.
[thetaB,rB,zB]=cart2pol(rbeadCutFormin(1,:),rbeadCutFormin(2,:),rbeadCutFormin(3,:));
xc1 = median(zB);
yc1 = r_ring-median(rB);
r = sqrt((zB-xc1).^2 + (r_ring-rB-yc1).^2);
dr = 0.01;
cnt = 1;
rList = (0:dr:2);
rCount = zeros(1,length(rList));
for k=rList
	cutCount = 0;
	for theta=linspace(0,2*pi,100)
		if yc1 + (k+dr/2)*sin(theta) < 0
			cutCount = cutCount + 1;
		end
	end
	rCount(cnt) = length(r(r >= k & r < k + dr))/((2*pi*(k+dr/2)*dr));
	cnt = cnt + 1;
end
normR = rCount/max(rCount);

% plot results, and find that r value where count drops to < 1% of maximum
% indR = find(normR < 0.03); %for delta myp2
indR = find(normR < 0.01);

% figure(1);
% plot(rList,normR,'bo');

% plot actin cross-section, with a black circle indicating where the bundle is.
% figure(2);
% plot(zB,r_ring-rB,'b.');
% hold on
% plot(xc1,yc1,'ko');
%rb = 0.095;
rb = rList(indR(1));
x = rb*sin(-pi:0.01*pi:pi) + xc1;
y = rb*cos(-pi:0.01*pi:pi) + yc1;
% c = [1 1 1];
% fill(x, y, c,'FaceColor','None','EdgeColor','k')
% axis equal;
% xlim([-0.3 0.3]);
% ylim([0 0.5]);
%xlim([-0.4 0.4]);
%ylim([0 0.7]);
% view(180,90); 

% calculate fraction of actin in bundle.
cnt = 0;
cntMyo = 0;
cntMyp = 0;
cntMypInBundle = 0;

[thetaT,rT,zT]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
[thetaT2,rT2,zT2]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));

rk = [];
zk = [];

for k=1:length(rB)

	if (zB(k) - xc1)^2 + (r_ring-rB(k)-yc1)^2 - rb^2 <= 0
		cnt = cnt + 1;
		
		rk = [rk r_ring-rB(k)];
		zk = [zk zB(k)];
		
		if (zB(k)-median(zT))^2 + (rB(k)-median(rT))^2 <= 0.051^2
			cntMyo = cntMyo + 1;
		end		
		
		if (zB(k)-median(zT2))^2 + (rB(k)-median(rT2))^2 <= 0.1^2
			cntMyp = cntMyp + 1;
		end
		
		%{
		for l=1:size(rmyo,2)
			[thetaT,rT,zT]=cart2pol(rmyo(1,l),rmyo(2,l),rmyo(3,l));
			if (zB(k)-zT)^2 + (rB(k)-rT)^2 <= 0.051^2
				cntMyo = cntMyo + 1;
				break;
			end
		end		
		
		for l=1:size(rmyp,2)
			[thetaT,rT,zT]=cart2pol(rmyp(1,l),rmyp(2,l),rmyp(3,l));
			if (zB(k)-zT)^2 + (rB(k)-rT)^2 <= 0.1^2
				cntMyp = cntMyp + 1;
				break;
			end
		end		
		%}
		
	end
end


[thetaB,rB,zB]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));

rk2 = [];
zk2 = [];
rmypBundle = [];

for k=1:size(rmyp,2)

	if (zB(k) - xc1)^2 + (r_ring-rB(k)-yc1)^2 - rb^2 <= 0
		cntMypInBundle = cntMypInBundle + 1;
		rk2 = [rk2 r_ring-rB(k)];
		zk2 = [zk2 zB(k)];
		rmypBundle = [rmypBundle rB(k)];
	end
	
end

sep = mean(rT) - mean(rmypBundle);

% disp('Statistics time')
% disp(['Fraction of actin in bundle % ' num2str(cnt/size(rbead,2)*100)])
% disp(['Fraction of myo2-actin in bundle % ' num2str(cntMyo/cnt*100)])
% disp(['Fraction of myp2-actin in bundle % ' num2str(cntMyp/cnt*100)])
% disp(['Fraction of myp2 in bundle % ' num2str(cntMypInBundle/size(rmyp,2)*100)])
% 
% % 95 thickness
% disp(['thickness ' num2str(prctile(rk,99) - prctile(rk,1))])
% disp(['width ' num2str(prctile(zk,99) - prctile(zk,1))])
% disp(['bundle size ' num2str(rb)])
% disp(['myp/myo separation ' num2str(sep)])

% 95 width

% fraction of beads in bundle int' myo2
