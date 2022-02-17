function [l nb] = getRandomL(vpol, rsev, kofffor, dbead)

%% Delta function distribution
l = 2.5 - vpol/kofffor/2;
% l = 2.4785;
% l = 1.7206;
nb = floor(l/dbead);
return;

%% Real distribution -- WARNING: CURRENTLY CHEATING TO FORBID SHORT FILS
% probability for a filament to have length lf
% consider filaments that are dbead to 100dbead long
%vpol = 0.1103;
%rsev = 0.0146;
%lfvec = 1:100;
%lfvec = lfvec * dbead;
%% probability 
%aa = (kofffor + lfvec * rsev)/vpol .* exp(-lfvec .* (2*kofffor + lfvec * rsev) / (2 * vpol));
%aa(1:nmin-1) = 0; % Forbid filaments with fewer than five beads
%aa = aa / sum(aa); % normalize the total to 1
%% cumulative probability
%acum = cumsum(aa);
%% number of actin beads not including formin
%nb = find(acum > rand, 1); 
%l = nb * dbead;
%return;
