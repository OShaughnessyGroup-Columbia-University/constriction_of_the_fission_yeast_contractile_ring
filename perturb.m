function vout = perturb(vin,head,eta)
% perturb SMOOTHLY every value of vin within a factor of eta, and keep the
% elements at positions "head" and "head - 1" fixed

%% generate a vector with npt elements, each with a value between (1-eta) and 1, and the first and last elements are 1
npt = 20;
pts = ones(1,npt);
for i = 2:npt-1
    pts(i) = 1 - rand * eta;
end
%% fit a smooth line to these points
x = linspace(0,1,npt);
f = fit(x', pts', 'smoothingspline','SmoothingParam', 1-1e-12);
%% apply the fitting to a vector of length vin
x = linspace(0,1,length(vin));
ff = f(x);
%% modify ff so that it starts and ends exactly at 1
% make the first point equal to 1
ff = ff + 1 - ff(1);
% scale linearly so that the end point is equal to 1
dif = 1 - ff(end);
for i = 2:length(ff) 
    ff(i) = ff(i) + dif * (i-1) / (length(ff)-1);
end
%% circularly shift the vector ff so that it starts at element # head
ff = circshift(ff,[0,head-1]);
%% multiply ff to vin
vout = vin .* ff';
%% shift vin back
vout = circshift(vout, [0,head]);
end