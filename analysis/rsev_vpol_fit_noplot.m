% close all
% disp("WARNING: INCREASING KOFFMYO TO ACCOUNT FOR LARGE EFFECTIVE FORMIN OFF RATE")
% koffmyo = koffmyo*1.75;
t_dissemble = 55; % 90% actin lost after 55 s treatment of lat A
lfvec = 0.1:0.1:10;

rsev_scan = linspace(1e-6,1.8/60*3,100); % avoid zero
vpol_scan = linspace(1e-6,0.24,100);

[rsev_grid,vpol_grid] = meshgrid(rsev_scan,vpol_scan);
l_mean = nan(size(rsev_grid));
ratio_t = nan(size(rsev_grid));

for i = 1:numel(rsev_grid)
    rsev = rsev_grid(i);
    vpol = vpol_grid(i);
    
    l_act_prob = l_act_dist(lfvec,rsev,vpol);
    l_mean(i) = sum(lfvec.*l_act_prob);
    l_t = l_act_t(lfvec,rsev,vpol,t_dissemble);
    ratio_t(i) = l_t / l_mean(i) * exp(-koffmyo * t_dissemble);
end
actin_per_fil = interp1([0,30],[500,97],real_t_min) * 350/200;% 500 at time 0, 97 at time 30 min
l_mean_goal = actin_per_fil * 500/190000;
ratio_t_goal = 0.1;
error2 = (1 - l_mean/l_mean_goal).^2 + (1-ratio_t/ratio_t_goal).^2;
[~,best] = min(error2(:));

% figure
% surf(rsev_grid,vpol_grid,l_mean)
% xlabel('r_{sev}')
% ylabel('v_{pol}')
% zlabel('mean l_{act} in ss.')
% hold on
% scatter3(1.8/60, .07, sum(lfvec.*l_act_dist(lfvec,1.8/60,0.07)),'ro')
% scatter3(rsev_grid(best),vpol_grid(best),l_mean(best),'yo')
% hold off

% figure
% surf(rsev_grid,vpol_grid,ratio_t)
% xlabel('r_{sev}')
% ylabel('v_{pol}')
% zlabel('l_{act}(t_dissemble) / l_{act, ss}.')
% hold on
% scatter3(1.8/60,.07,0.1,'ro')
% scatter3(rsev_grid(best),vpol_grid(best),ratio_t(best),'yo')
% hold off

rsev = rsev_grid(best);
vpol = vpol_grid(best);

if(rsev == 1e-6)
    vpol = l_mean_goal*koffmyo;
end

% disp("WARNING: RESTORING KOFFMYO")
% koffmyo = koffmyo/1.75;