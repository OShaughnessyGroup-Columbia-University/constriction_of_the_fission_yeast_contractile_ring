[theta2,r2,z2]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));
[theta4,r4,z4]=cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));
[theta3,r3,z3,flg]=cutWhiskers(theta4,r4,z4,r_ring);
naiveProjectedReachFlag = zeros(1,length(r3));
reachFlag = zeros(1,length(r3));
rcut = 0.1;

for k=1:length(r2)
	for l=1:length(r3)
		if ((r2(k)-r3(l)).^2 + (z2(k)-z3(l)).^2 <= rcut.^2)
			naiveProjectedReachFlag(l) = naiveProjectedReachFlag(l) + 1;
		end
		if ((rmyp(1,k)-rbead(1,l)).^2 + (rmyp(2,k)-rbead(2,l)).^2  + (rmyp(3,k)-rbead(3,l)).^2  <= rcut.^2)
			reachFlag(l) = reachFlag(l) + 1;
		end
	end
end

%{
plot3(rbead(1,reachFlag > 0 & flg),rbead(2,reachFlag > 0 & flg),rbead(3,reachFlag > 0 & flg),'b.')
hold on
plot3(rbead(1,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),rbead(2,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),rbead(3,naiveProjectedReachFlag > 0 & reachFlag == 0 & flg),'g.')
hold on
plot3(rmyp(1,:),rmyp(2,:),rmyp(3,:),'r.')
axis equal;
%}