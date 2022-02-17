clear
t = 1:26;
time_wu = [-5, 5, 15, 25];
rho_myo2_exp = [18, 22, 32, 52];
rho_myp2_exp = [9, 17, 32, 34];
rhom = spline(time_wu, rho_myo2_exp, t)/18 * 290 + spline(time_wu, rho_myp2_exp, t)/9 * 200;
l = interp1([0,30],[500,97],t)*350/200;

koffmyo = 0.0245;

[~,vpol] = rsev_vpol_fit_noplot_fn(l,koffmyo);

vmyo0 = 0.24;
fstall = 1;
tens = .5 * fstall * rhom .* (l* 500/190000) .* (1- vpol/vmyo0); 

phi = acosd((1.85-.07*t)/1.85);
plot(phi,tens)
axis([-Inf,Inf,0,Inf])