% Note: this script is still in development.
% please be careful while using it.

% all units are pN,s,micron unless otherwise stated.
% all positions are relative to myo2 center except
% when results are plotted, where y = 0 is the septum 
% i.e. the wall.

% set a seed
rng(10);

Nbeads = 40;					% number of beads in cross-section
rbeads = randn(2,Nbeads).*0.05; % x,y positions of beads in c.s.

dt = 0.1; 			% seconds, timestep of simulation

gammaAct = 300; 	% pN s/micron, drag coefficient of motion in cross-section
					% in the simulation. per bead drag coefficient is about
					% 20, but it is attached to a filament, so the gammaAct
					% used here has to be somewhat larger. It is unclear
					% what exact value I should use here.

gammaMyp = 20;		% pN s/micron, drag coefficient of myp2 in cross-section
					% note: as myp2 motion is not properly implemented, I don't
					% know how large this is.
					
T = 1000; 			% total time of simulation. Strictly speaking, this
					% should be chosen such that the velocities are very small.

rexc = 0.02/sqrt(2);% excluded volume radius of actin
f0exc = 20; 		% maximum actin excluded volume force

ymyp = -27/1000;	% y position of myp2
yfor = 50/1000;		% y position of formin, which is fixed.
flaplace = 0.6;		% inward laplace force per rod (100 nm of actin)

% force constant of wall.
kWall = 3000;
rcutWall = 100/1000;
yWall = +94/1000;	% position of wall

% node-myp2 excluded volume parameters
% measure from center of the formin
rexPO = 202/1000;
kexPO = 70;

% plot ring every
plotEvery = 200;

cnt = 0;

% perform simulation
for k=1:round(T/dt)
	[vact,vmyp] = getVelocity(rbeads,ymyp,gammaAct,gammaMyp,...
								rexc,f0exc,flaplace,kWall,rcutWall,...
								rexPO,kexPO,yfor,yWall);
	rbeads = rbeads + vact*dt;
	
	% Myp2 motion part of the code under development.
	% ymyp = ymyp + vmyp*dt;
    cnt = cnt + 1;
    
    if cnt > plotEvery
        cnt = 0;
        clf;

        % IMPORTANT: everything is calculated assuming myo2
        % is fixed at x,y = 0,0. Now we shift everything downward

        % plot results
        viscircles([rbeads(1,:);rbeads(2,:)-yWall]', rexc*ones(1,Nbeads),'Color','k');
        hold on;

        % plot saucepans of myosins
        plot(0.025*cos([0:100]*2*pi/100),0.025*sin([0:100]*2*pi/100)-yWall,'m-');
        plot(0.05*cos([0:100]*2*pi/100),0.05*sin([0:100]*2*pi/100)-yWall,'r-');

        plot(0.05*cos([0:100]*2*pi/100),0.05*sin([0:100]*2*pi/100)+ymyp-yWall,'c-');
        plot(0.1*cos([0:100]*2*pi/100),0.1*sin([0:100]*2*pi/100)+ymyp-yWall,'b-');

        % wall position
        hline(0)
        xlim([-0.2 0.2]);
        ylim([-0.2 0.2]);
        axis equal;
        
        pause(0.5);
    end
    
end

figure;
histogram(sqrt(sum(vact.*vact,1))*1000);
xlabel('Actin bead velocities (nm/s)')
ylabel('Count');

function [vact,vmyp] = getVelocity(x,ymyp,gammaAct,gammaMyp,...
									rexc,f0exc,flaplace,...
									kWall,rcutWall,rexPO,kexPO,...
									yfor,yWall)

	N = size(x,2);

	% apply myo2 saucepan
	r = sqrt(sum(x.^2,1));
	r(r < 1e-3) = 1e-3;
	rhat = x./r;
	fo = (fmyo2(r)).*(rhat);
	
	% apply myp2 saucepan
	y = x;
	y(2,:) = y(2,:) - ymyp;
	r = sqrt(sum(y.^2,1));
	r(r < 1e-3) = 1e-3;
	rhat = y./r;
	fp = (fmyp2(r)).*(rhat);
	
	% apply laplace force;
	yhat = rhat*0;
	yhat(2,:) = -1;
	fl = flaplace * yhat;
	
	% apply excluded volume forces
	fe = fl * 0;
	sigma_sq = rexc^2;
	
	% loop over all pairs of beads
	for k=1:N
		for l=1:N
		
			% if same bead, continue
			if k == l
				continue
			end
			
			% vector between beads and distance
			nhat = [x(1,k)-x(1,l);x(2,k)-x(2,l)];
			dst = norm(nhat);
			nhat = nhat./dst;
			
			% calculate forces
			fe(:,k) = fe(:,k) + (f0exc * exp(-dst.^2/(2*sigma_sq))) * nhat;
		end
	end
	
	% calc velocity
	vact = (fo + fp + fl + fe)./gammaAct;
	
	%{
	% Myp2 motion part of the code under development.
	
	% myp2 excluded volume force
	fexPO = 0;
	
	if((abs(ymyp-yfor)) < rexPO )
		fexPO = kexPO * (rexPO - abs(ymyp-yfor));
		if ymyp < yfor
			fexPO = -fexPO;
		end
	end
	
	% calculate myp2 velocity
	fcapReaction = -sum(fp,2);
	vmyp = (fcapReaction(2)-kWall*abs(ymyp-yWall+rcutWall)+fexPO)./gammaMyp;
	%}
	vmyp = 0;
	
end

function b = fmyo2(r)
	% myo2 saucepan force
	b = 0*r;
	b(r >= 0.025) = -7-7/(0.025) * (r(r >= 0.025) - 0.025);
	b(r >= 0.05) = 0;
end

function b = fmyp2(r)
	% myp2 saucepan force
	b = 0*r;
	b(r >= 0.05) = -7-7/(0.05) * (r(r >= 0.05) - 0.05);
	b(r >= 0.1) = 0;
end