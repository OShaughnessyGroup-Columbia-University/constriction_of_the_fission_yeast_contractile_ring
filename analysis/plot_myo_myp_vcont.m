% function plot_myo_myp_concentration
figure
tlist = 0:10:1200;
clear nmyo nmyp nact
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_1_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_1_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('matur')

figure
clear nmyo nmyp nact
tlist = 0:10:1200;
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_2_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_2_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('assmyp')

figure
clear nmyo nmyp nact
tlist = 0:10:1200;
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_3_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_3_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('allmyo')

figure
clear nmyo nmyp nact
tlist = 0:10:660;
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_4_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_4_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('rass')

figure
clear nmyo nmyp nact
tlist = 0:10:1200;
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_5_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_5_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('allmyo old')

figure
clear nmyo nmyp nact
tlist = 0:10:1200;
for i = 1:length(tlist)
    if(exist(['vcont2/vcont2_6_' num2str(tlist(i)) 'sec.mat']))
        load(['vcont2/vcont2_6_' num2str(tlist(i)) 'sec.mat']);
    end
    nmyo(i) = length(rmyo);
    nmyp(i) = length(rmyp);
    nact(i) = length(rbead);
end
plot(tlist(1:length(nmyo)), nmyo, tlist(1:length(nmyo)), nmyp, tlist(1:length(nmyo)), nact/100)
ylim([0, 250])
title('allmyo old')