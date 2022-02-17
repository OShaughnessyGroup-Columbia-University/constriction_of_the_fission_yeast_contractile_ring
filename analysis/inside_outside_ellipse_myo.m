[theta2,r2,z2]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
rb = rbead;
rb(:,ifor) = [];
bnd = logical(bundle(rb));
[theta3,r3,z3]=cart2pol(rb(1,bnd),rb(2,bnd),rb(3,bnd));
% [theta3,r3,z3,flg]=cutWhiskers(theta4,r4,z4,r_ring);
% length(find(flg))/length(flg)
theta3 = mod(theta3, 2*pi);
naiveProjectedReachFlag = zeros(1,length(r3));
reachFlag = zeros(1,length(r3));
rcut_b = 0.051;
rcut_a = 0.066;

for k=1:length(r2)
	for l=1:length(r3)
		if ((r2(k)-r3(l)).^2 + (z2(k)-z3(l)).^2 <= rcut_b.^2)
			naiveProjectedReachFlag(l) = naiveProjectedReachFlag(l) + 1;
		end
		r = rmyo(:,k)-ract(:,l);
		thetaHat = [-rmyo(2,k),rmyo(1,k),0];
		thetaHat = thetaHat./sqrt(sum(thetaHat.^2));
		if ((dot(r,thetaHat))^2 * (1-rcut_b^2/rcut_a^2) >= sum(r.^2) - rcut_b.^2)
			reachFlag(l) = reachFlag(l) + 1;
            break
		end
	end
end
% scatter( phiact(reachFlag==0), z3(reachFlag==0))
% scatter( phiact(reachFlag>0), z3(reachFlag>0))

% scatter( phiact(reachFlag==0), r3(reachFlag==0))
scatter( phiact(reachFlag>0),  r3(reachFlag>0))

%{
plot3(rbead(1,(reachFlag > 0) & flg),rbead(2,(reachFlag > 0) & flg),rbead(3,(reachFlag > 0) & flg),'b.')
hold on
plot3(rbead(1,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),rbead(2,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),rbead(3,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),'g.')
hold on
plot3(rmyo(1,:),rmyo(2,:),rmyo(3,:),'r.');
axis equal;
%}