% function plot_myo_myp_concentration
figure
tlist = 0:10:1200;
clear nmyo nmyp nact nfor
for i = 1:length(tlist)
    load(['matur2/matur2_1_' num2str(tlist(i)) 'sec.mat']);
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
    nfor(i) = length(ifor);
    if r_ring < 1.85
        break;
    end
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)),...
    nact/100, tlist(1:length(nmyo)), nfor)
ylim([0, 250])
title('matur')
legend({'myo', 'myp', 'actin', 'formin'})

figure
clear nmyo nmyp nact nfor
tlist = 0:10:1200;
for i = 1:length(tlist)
    load(['mypdelay/mypdelay_1_' num2str(tlist(i)) 'sec.mat']);
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
    nfor(i) = length(ifor);
    if r_ring < 1.85
        break;
    end
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)),...
    nact/100, tlist(1:length(nmyo)), nfor)
ylim([0, 250])
title('mypdelay')
legend({'myo', 'myp', 'actin', 'formin'})

figure
clear nmyo nmyp nact nfor
tlist = 0:10:2400;
for i = 1:length(tlist)
    load(['allmyo3/allmyo3_1_' num2str(tlist(i)) 'sec.mat']);
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
    nfor(i) = length(ifor);
    if r_ring < 1.85
        break;
    end
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)),...
    nact/20, tlist(1:length(nmyo)), nfor)
ylim([0, 250])
title('allmyo')
legend({'myo', 'myp', 'actin', 'formin'})

figure
clear nmyo nmyp nact nfor
tlist = 0:1:22;
for i = 1:length(tlist)
    load(['wtfc4/wtfc4_1_' num2str(tlist(i)) 'min.mat']);
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
    nfor(i) = length(ifor);
%     if r_ring < 1.85
%         break;
%     end
end
tlist = 60*tlist;
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)),...
    nact/20, tlist(1:length(nmyo)), nfor)
ylim([0, 250])
title('wtfc3')
legend({'myo', 'myp', 'actin', 'formin'})

% figure
% clear nmyo nmyp nact
% tlist = 0:10:660;
% for i = 1:length(tlist)
%     load(['rass/rass_2_' num2str(tlist(i)) 'sec.mat']);
%     nmyo(i) = length(rmyo);
%     nmyp(i) = length(rmyp);
%     nact(i) = length(rbead);
%     if r_ring < 1.85
%         break;
%     end
% end
% plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
% ylim([0, 250])
% title('rass')
% 
% figure
% clear nmyo nmyp nact
% tlist = 0:10:1200;
% for i = 1:length(tlist)
%     load(['allmyo_old/allmyo_1_' num2str(tlist(i)) 'sec.mat']);
%     nmyo(i) = length(rmyo);
%     nmyp(i) = length(rmyp);
%     nact(i) = length(rbead);
%     if r_ring < 1.85
%         break;
%     end
% end
% plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
% ylim([0, 250])
% title('allmyo old')
% 
% figure
% plot(tlist(1:length(nmyo)), nmyp./nact/100)