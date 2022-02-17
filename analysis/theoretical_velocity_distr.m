%% set up the pdf of filament lengths
%     the functional form of the pdf below roughly
%     matches the distribution that zach sees
%     in the simulation.

kofffor = 1/40.;
% rsev = 0.0155;
% vpol = 0.127;
gamma_anc = 1500;
theoretical_filament_histogram
% for k=1:length(l)
%     p(k)=abs(exp(-(l(k)-1).^2 ./7) - exp(-l(k)./0.2));
%     
% end
l = 0.01:0.01:10;
p = histcounts(Lflts, [0 l]);
% p = (kofffor + l * rsev)/vpol .* exp(-l .* (2*kofffor + l * rsev) / (2 * vpol));
% p(l == 1.0) = 1;
% p(l < 0.01) = 0; % get rid of filaments shorter than lmin um

%{
for k=1:length(l)
	if (l(k) < 3)
		p(k)=1;
	else
		p(k)=exp(-(l(k)-3));
	end
	
	if(l(k) < 0.6)
		p(k) = l(k)/0.6;
	end
end
%}
p = p/sum(p);
% sum(p.*l)

% figure
% plot(l, p*0.01)
% xlabel('Filament length (\mum)')
% ylabel('Probability density')
% set(gca, 'Linewidth', 3)
% set(gca, 'FontSize', 24)
% set(gca, 'Position', [0.2 0.2 0.7 0.7]);

%% set up params
%    not sure the exact vals of params below matters -
%    I think it just rescales the velocity, and the mean
%    of the velocity can be fit to any value by tuning these
%    parameters.

f = 2.*(1 - vpol/0.24); % force per head of myosin
rho_head = (2900+2000)/500; % density of heads per micron 
		                    % of filament (assume Myo2 = Myp2)
% gamma_anc = 1000;      % node anchor drag coefficient

% number of filaments in the core
Ncore = 30;

% set up arrays
N = 1e6;
v = zeros(1,N);
ldist = zeros(1,N);

% sample the l distribution
cnt = 1;

for k=1:length(l)
    M = round(p(k)*N);
    ldist(cnt:min(cnt + M,N)) = l(k);
    cnt = cnt + M + 1;
end

ldist=ldist(randperm(length(ldist)));
pdist = sign(rand(1, length(ldist))*2-1);

%% get the number of filaments per node
% all nodes have 1 formin
% Nfildist = ones(1,N);


% basically, every incoming formin has a 
% probability to bind a node. this probability
% is proportional to (4-n) where n is the number
% of formins already in the node.

Nfildist = zeros(1,N);
Nfor = round(200/180 * N);
pArray = [4 3 2 1 0]./10;
%pArray = [8 2 0 0 0]./10;

for k=1:Nfor
	success = false;
	while ~success
		cnt = randi(length(Nfildist));
		Nfilcurr = Nfildist(cnt);
		pChoose = pArray(Nfilcurr + 1);
		success = (rand() < pChoose);
	end
	Nfildist(cnt) = Nfildist(cnt) + 1;
end

%% get the velocities
cnt = 1;

for k=1:N

    x = round(rand(Ncore,1)); % ~Ncore filaments in the core
    x = x*2-1;   % + or -1 polarity
    P = sum(x);  % net polarity
	
	Nfil = Nfildist(k);
	
    v(k) = (P*16/Ncore); 
		% 16/Ncore multiplies the polarity as there are 16 heads
		% per cluster and ~Ncore flts in the core.
		
	for il = 1:Nfil
		v(k) = v(k) + pdist(cnt)*rho_head*ldist(cnt);
		cnt = cnt + 1;
		if cnt > N
			break;
		end
	end
		
	v(k) = v(k)*(f/gamma_anc); 
	
	if cnt > N
		break;
	end
end

%% histogram of node velocities excluding nodes with no formins
figure
hold on
v = v(1:k);
Nfildist = Nfildist(1:k);
% account for 20% of whiskers and other noise
%v = v + randn(1, length(v))*10/gamma_anc;
% bb = rand(1,length(v)) > 0.8;
% v(bb) = randn(1,sum(bb))*10/gamma_anc;
for nfmin = 1:4
%     v = v(Nfildist > nfmin);
    vmax = max(abs(v))*1e3;
    histogram(v(Nfildist >= nfmin)*1000, -vmax:4:vmax, 'Normalization', 'pdf', 'Displayname', ['nf>=' num2str(nfmin)]);
end
xlabel('Formin velocity (nm/s)')
ylabel('Probability density')
title(['Theoretical v_{for}'])
set(gca, 'Linewidth', 3)
set(gca, 'FontSize', 24)
set(gca, 'Position', [0.2 0.2 0.7 0.7]);
cnts = histcounts(v*1000, -vmax:2:vmax);
peak = max(cnts);
% p20_ratio = peak/cnts(floor(length(cnts)/2))

% %% formin velocity distribution
% figure
% hold on
% box on
% pldist = pdist.*ldist;
% for nfor = 1:4
%     nval = floor(length(ldist)/nfor);
%     plmat = [];
%     for i = 1:nfor
%         plmat = [plmat; pldist(1+(i-1)*nval:i*nval) ];
%     end
%     histogram(sum(plmat, 1), -15:0.25:15, 'Normalization', 'probability')
% end
% xlabel('Net length of actin (\mum)')
% ylabel('Probability')
% set(gca, 'Linewidth', 3)
% set(gca, 'FontSize', 24)
% set(gca, 'Position', [0.2 0.2 0.7 0.7]);