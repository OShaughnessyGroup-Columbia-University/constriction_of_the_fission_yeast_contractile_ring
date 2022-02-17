% clear
% load shrinkp_48_0min
% load shrinkp_48_2min
clear
irun = 4;
figure
load(['gm_2/gm_2_' num2str(irun) '_00sec.mat'])
load(['gm_2/gm_2_' num2str(irun) '_600sec.mat'])
% nstep = 1;
% dt = 0.2; % set to 0.2 to avoid update
% rk1_v32
% disp(size(vmyo))
% myo2;
phimyo = cart2pol(rmyo(1,:), rmyo(2,:));
phihat = [-sin(phimyo) ; cos(phimyo); zeros(1,length(phimyo))];
for im=1:length(vmyo)
    vphi(im) = phihat(:,im)'*vmyo(:,im);
end

myo2_speed = sqrt(sum(vmyo.*vmyo));
myo2_sp_mean=mean(myo2_speed);
edges = -.1:.002:0.1;
h1=histogram(vphi,edges);
% mean(abs(vphi))
% [n,edges]=histcounts(myo2_speed);
% x = (edges(1:end-1)+edges(2:end))/2;
% stem(x,n,'r-')
hold on
% formin
% formin_speed = sqrt(sum(vbead(:,ifor).*vbead(:,ifor)));
% h2=histogram(formin_speed,'FaceColor','blue');
% % myp2
% % myp2_speed = sqrt(sum(vmyp.*vmyp));
% % h3=histogram(myp2_speed,'FaceColor','green');
% % [n,edges]=histcounts(myp2_speed);
% % x = (edges(1:end-1)+edges(2:end))/2;
% % stem(x,n,'g-')
% h1.Normalization = 'probability';
% % h1.BinCounts = 30;
% h2.Normalization = 'probability';
% h2.BinWidth = h1.BinWidth;
% h3.Normalization = 'probability';
% h3.BinWidth = 0.05;
% hold off
xlabel('Speed (\mum/s)')
ylabel('Frequency')
% legend('Myo2p','Formin','Location','northeast')