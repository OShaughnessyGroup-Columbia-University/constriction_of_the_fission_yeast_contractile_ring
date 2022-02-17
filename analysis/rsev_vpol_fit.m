close all
kforoff = 0.0245;
t = 55; % 90% actin lost after 55 s treatment of lat A
lfvec = 0.1:0.1:10;

rsev_scan = linspace(1e-5,1.8/60*3,100); % avoid zero
vpol_scan = linspace(1e-5,0.07*3,100);

[rsev_grid,vpol_grid] = meshgrid(rsev_scan,vpol_scan);
l_mean = nan(size(rsev_grid));
ratio_t = nan(size(rsev_grid));

for i = 1:numel(rsev_grid)
    rsev = rsev_grid(i);
    vpol = vpol_grid(i);
    
    l_act_prob = l_act_dist(lfvec,rsev,vpol);
    l_mean(i) = sum(lfvec.*l_act_prob);
    l_t = l_act_t(lfvec,rsev,vpol,t);
    ratio_t(i) = l_t / l_mean(i) * exp(-kforoff * t);
end

l_mean_goal = 2.5;
ratio_t_goal = 0.1;
error2 = (1 - l_mean/l_mean_goal).^2 + (1-ratio_t/ratio_t_goal).^2;
[~,best] = min(error2(:));

surf(rsev_grid,vpol_grid,l_mean)
xlabel('r_{sev}')
ylabel('v_{pol}')
zlabel('mean l_{act} in ss.')
hold on
scatter3(1.8/60,.07,sum(lfvec.*l_act_dist(lfvec,1.8/60,0.07)),'ro')
scatter3(rsev_grid(best),vpol_grid(best),l_mean(best),'yo')
hold off

figure
surf(rsev_grid,vpol_grid,ratio_t)
xlabel('r_{sev}')
ylabel('v_{pol}')
zlabel('l_{act}(t) / l_{act, ss}.')
hold on
scatter3(1.8/60,.07,0.1,'ro')
scatter3(rsev_grid(best),vpol_grid(best),ratio_t(best),'yo')
hold off

rsev_grid(best)
vpol_grid(best)