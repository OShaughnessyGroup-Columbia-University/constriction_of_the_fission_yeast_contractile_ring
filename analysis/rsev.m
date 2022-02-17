t = 60;
kofffor = 0.0052; 
vp = .07;
rsev = .05;
fn = @(x) x/(1+.5*rsev);
integral(fn, 0, 1)