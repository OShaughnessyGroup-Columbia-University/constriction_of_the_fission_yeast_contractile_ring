% clear
% load shrinkp_48_0min
% load shrinkp_48_2min
nstep = 1;
dt = 0.2; % set to 0.2 to avoid update
rk1_v12
% disp(size(vmyo))
% myo2
myo2_speed = sqrt(sum(vmyo.*vmyo));
h1=histogram(myo2_speed,30,'FaceColor','red');
% [n,edges]=histcounts(myo2_speed);
% x = (edges(1:end-1)+edges(2:end))/2;
% stem(x,n,'r-')
hold on
% formin
formin_speed = sqrt(sum(vbead(:,ifor).*vbead(:,ifor)));
h2=histogram(formin_speed,'FaceColor','blue');
% myp2
% myp2_speed = sqrt(sum(vmyp.*vmyp));
% h3=histogram(myp2_speed,'FaceColor','green');
% [n,edges]=histcounts(myp2_speed);
% x = (edges(1:end-1)+edges(2:end))/2;
% stem(x,n,'g-')
h1.Normalization = 'probability';
% h1.BinCounts = 30;
h2.Normalization = 'probability';
h2.BinWidth = h1.BinWidth;
% h3.Normalization = 'probability';
% h3.BinWidth = 0.05;
hold off
xlabel('Speed (\mum/s)')
ylabel('Frequency')
legend('Myo2p','Formin','Location','northeast')