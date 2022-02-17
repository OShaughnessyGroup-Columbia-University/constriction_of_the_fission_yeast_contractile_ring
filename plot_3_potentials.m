temp = linspace(-0,200,200);

% x: distance between myo2 and myp2
% y: distance between myo2 and actin
[x,y] = meshgrid(temp,temp);

rcap_myo = 51;
rcap_myp = 100;

% between myo2 and myp2
prefactor = .3;
kexc = prefactor * 15 / 1000 * 178*200;
u1 = .5 * kexc * (x - rcap_myo - rcap_myp) .^ 2;
u1(x - rcap_myo - rcap_myp > 0) = 0;

% between myo2 and actin
k_myo = 100 / 1000 * 178 * 50;
u2 = .5 * k_myo * y.^2;
u2(abs(y) > rcap_myo) = .5 * k_myo * rcap_myo^2;

% between myp2 and actin
k_myp = 25 / 1000 * 200 * 50;
u3 = .5 * k_myp * (x-y) .^ 2;
u3(abs(x-y) > rcap_myp) = .5 * k_myp * rcap_myp^2;

% equivalent potential energy due to ring tension
tens = 0;
u4 = - tens * 2 * pi * y;

% total potential energy
utot = u1+u2+u3+u4;

% find minimum of energy
[c i] = min(utot(:));

plot3(x(i),y(i),utot(i),'ro')
hold on
surf(x,y,utot,'EdgeColor','none')
hold off
xlabel('myo2-myp2 distance (nm)')
ylabel('myo2-actin distance (nm)')
zlabel('Potential (A.U.)')
view(0,90)

% limit on kexc
alpha = (rcap_myo / (rcap_myo+rcap_myp))^2;
first_limit = alpha * k_myo * k_myp / ((1-alpha) * k_myp - alpha * k_myo)

beta = (rcap_myp / (rcap_myo+rcap_myp))^2;
second_limit = beta * k_myo * k_myp / ((1-beta) * k_myo - beta*k_myp)

kexc_limit = min(first_limit, second_limit)

% limit on x (myo2-myp2 distance)
max_myo_myp_dist = kexc_limit*(k_myo+k_myp)/(k_myo*k_myp + kexc_limit*(k_myo+k_myp)) * (rcap_myo+rcap_myp)