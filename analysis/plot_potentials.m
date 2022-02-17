x = linspace(-200,200,200);
c_myp = 0;
c_myo = 70;
rcap_myo = 51;
rcap_myp = 70;
k_myo = 2;
k_myp = 1;
u_myo = .5 * k_myo * (x-c_myo).^2;
u_myo(abs(x-c_myo)>rcap_myo) = .5 * k_myo * rcap_myo^2;
plot(x,u_myo,'b-')
hold on
u_myp = .5 * k_myp * (x-c_myp).^2;
u_myp(abs(x-c_myp)>rcap_myp) = .5 * k_myp * rcap_myp ^2;
plot(x,u_myp,'r-')
plot(x,u_myo+u_myp,'k-','LineWidth',2)
hold off
legend('myo2','myp2','total')
xlabel('x (nm)')
ylabel('Potential (A.U.)')