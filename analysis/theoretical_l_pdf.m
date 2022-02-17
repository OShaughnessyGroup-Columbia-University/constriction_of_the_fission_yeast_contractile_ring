% Script to calculate the theoretical pdf of filament lengths
% and some associated quantities.

% i.e. calculate p.d.f. of filament lengths over windows of size
% 0,T,2T etc. where T is a timescale. Over a time of zero, the pdf
% will just be the pdf calculated by Joseph (& others, this is an old result).
% Over longer intervals, the pdf can be calculated by evolving the 'joseph' pdf
% over time, switching off node turnover.

% parameters of the model: node turnover rate, actin growth and severing rate.
% base units are micron,second,picoNewton.

kofffor = 1/42;
vpol = 127/1000;
rsev = 0.0155; % this is per micron per second.

% calculate pdfs with a bin size of 100 nm.
dl = 1/10;
l = 0:dl:10;

% this is the ideal pdf of instantaneous filament lengths
% i.e. the 'joseph' pdf
%plot(l,idealPdf(l,kofffor,vpol,rsev),'b-');
%hold on;

% this is the ideal pdf of a window of size infinity.
plot(l,idealPdf(l,0*kofffor,vpol,rsev),'Color','#77AC30');
hold on;

% generate filament lengths from the 'joseph' pdf
l = rand(1,1000000) * 10;
lList = l(rand(1,length(l)) < idealPdf(l,kofffor,vpol,rsev));
histogram(lList,'Normalization','pdf','DisplayStyle','Stairs','EdgeColor','red')
mean(lList)
lList0 = lList;

% evolve them over different windows, and plot the mean
% filament length as measured over the respective windows.

dt = 0.01;
Ttot = 10;
cnt = 1;

for k=1:6
	for t=0:dt:Ttot
		lList = lList + vpol*dt;
		toBeCut = (rand(1,length(lList)) < dt*rsev*lList);
		lList(toBeCut) = rand(1,length(lList(toBeCut))) .* lList(toBeCut);
		lList0 = lList0 + lList;
		cnt = cnt + 1;
	end
	mean(lList)
	histogram(lList0/cnt,'Normalization','pdf','DisplayStyle','Stairs','EdgeColor','black');
end

legend('ideal inf','0','10','20','30','40','50','60');

function p = idealPdf(l,kofffor,vpol,rsev)
	% calculate the 'joseph' pdf i.e. steady-state
	% pdf of a bunch of filaments undergoing growth,
	% severing, and turnover.
	p = ((kofffor + l * rsev)/vpol) .* exp(-l .* (2*kofffor + l * rsev) ./ (2 * vpol));
end