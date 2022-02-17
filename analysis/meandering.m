load fh2d9_42_00min
load fh2d9_42_9min

[~,rho] = cart2pol(rbead(1,:),rbead(2,:));
d = ((1.85-0.07*9) - rho)*1000;

subplot(2,2,1)
hold on
for i = randsample(length(ifor),10)'
    plot(d(ifor(i):ipt(i)))
end
xlabel('bead index from barbed end')
ylabel('distance from membrane (nm)')

plot([0,50],[141,141],'k--')

[~,temp] = cart2pol(rmyp(1,:),rmyp(2,:));
mean_myp2_d = ((1.85-0.07*9) - mean(temp))*1000;

plot([0,50],[mean_myp2_d+100,mean_myp2_d+100],'k--')

hold off
axis([0,30,0,Inf])

subplot(2,2,2)
hold on
for i = randsample(length(ifor),10)'
    plot(autocorr(d(ifor(i):ipt(i))))
end
xlabel('bead index')
ylabel('ACF')
hold off

subplot(2,2,3)
c = zeros(size(d));
c(d <= 141) = -1;   % first class, bound to both Myo2 and Myp2
c(d > 141) =  1;   % second class, not bound to Myo2
hold on
for i = randsample(length(ifor),10)'
    plot(autocorr(c(ifor(i):ipt(i))))
end
xlabel('bead index')
ylabel('ACF')
hold off

subplot(2,2,4)
acf = [];
mu = mean(c);
for i = 0:30
    x = nan(1,size(rbead,2));
    x_ind = 1;
    for j = 1:length(ifor)
        for k = ifor(j):(ipt(j)-i)
            x(x_ind) = (c(k)-mu)*(c(k+i)-mu);
            x_ind = x_ind + 1;
        end
    end
    exx = nanmean(x);
    acf = [acf, exx/var(c)];
end
plot(0:30,acf)
xlabel('bead index')
ylabel('ACF')
