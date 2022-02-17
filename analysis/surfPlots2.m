[theta1,r1,z1]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));
[theta2,r2,z2]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
[theta3,r3,z3]=cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));

z = [z1,z2,z3];
M = zeros(1000*length(z),2);
t = [ones(1,length(z1)),2*ones(1,length(z2)),3*ones(1,length(z3))];
wdth = [0.1*ones(1,length(z1)),0.05*ones(1,length(z2)),0.0025*ones(1,length(z3))];

for k=1:length(z)
	M(1000*(k-1)+1:1000*(k-1)+1000,1) = random('Uniform',z(k)-wdth(k),z(k)+wdth(k),[1000,1]);
	M(1000*(k-1)+1:1000*(k-1)+1000,2) = t(k);
end

%{
M = zeros(length(z1) + length(z2) + length(z3),2);

t1 = 1:length(z1);
M(t1,1) = z1;
M(t1,2) = 0;

t1 = length(z1)+1:length(z1)+length(z2);
M(t1,1) = z2;
M(t1,2) = 1;

t1 = length(z1)+length(z2)+1:length(z1)+length(z2)+length(z3);
M(t1,1) = z3;
M(t1,2) = 2;
%}