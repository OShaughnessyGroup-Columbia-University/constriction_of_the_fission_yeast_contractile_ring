phi_bleb = [20.49229194;
    25.9006289423;
    21.1563478959;
    3.8384221977;
    29.0359422986;
    17.3831358875;
    7.9229124337;
    33.3860885112;
    19.2633525967]


t_bleb = [524.1307652329;
    476.646624327;
    373.5907863905;
    78.42257185;
    466.7396478697;
    178.6180265754;
    189.2669551875;
    417.4876347527;
    301.2296062631];

t_bleb = [phi_bleb, t_bleb];

figure
phi = 10:10:80;
for i = 1:8
    n_bleb_bin(i) = length(find(discretize(t_bleb(:, 1) , 5:10:65)==i));
    t_bleb_bin(i) = mean(t_bleb(discretize(t_bleb(:, 1) ,...
        5:10:65)==i , 2));
    t_bleb_bin_std(i) = std(t_bleb(discretize(t_bleb(:, 1) , ...
        5:10:65)==i , 2));
end
errorbar(phi, t_bleb_bin, t_bleb_bin_std, 'o', 'MarkerSize', 8, 'LineWidth', 2)