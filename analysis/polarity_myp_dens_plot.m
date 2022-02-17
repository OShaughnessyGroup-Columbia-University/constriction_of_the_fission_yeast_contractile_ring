figure
hold on
%% myp histogram
phi_bin = -pi:pi/10:pi;
[phimyp rhomyp zmyp] = cart2pol(rmyp(1, :), rmyp(2, :), rmyp(3, :));
myp_dens = histcounts(phimyp, phi_bin);
plot((phi_bin(1:end-1) + phi_bin(2:end))/2, myp_dens)
ylim([0, 20])

%% getting polarities
rdiff = get_rdiff(rbead, 1, ipt);
phibead = cart2pol(rbead(1, :), rbead(2, :), rbead(3, :));
phihat = [-sin(phibead); cos(phibead); zeros(1, length(phibead))];
for i = 1:length(phibead)
    p(i) = sign(sum(phihat(:, i).*rdiff(:, i)));
end

%% binning
yyaxis right
ploc = zeros(1, length(phi_bin)-1);
for jbead = 1:length(phibead)
    for i=1:length(phi_bin)-1
        if(phi_bin(i) < phibead(jbead) && phibead(jbead) < phi_bin(i+1))
%             disp(phi_bin(i))
%             disp(phibead(jbead))
%             disp (phi_bin(i+1))
            ploc(i) = ploc(i) + p(jbead);
        end
    end
end
plot((phi_bin(1:end-1) + phi_bin(2:end))/2, ploc)
