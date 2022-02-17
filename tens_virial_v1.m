%% Calculate tension by centripetial forces on the anchors
if mod(length(fanc),3) ~=0
    error('the length of fc is not a multiple of 3')
end
fc3 = vec2mat(fc,3)';
if mod(length(force),3) ~=0
    error('the length of force is not a multiple of 3')
end
force_vec = vec2mat(force,3)';
f_tot = fc3 + force_vec;
r_all = [rbead,rmyo,rmyp];
[th_all,~,rho_all] = cart2sph(r_all(1,:), r_all(2,:), r_all(3,:));

r_all(3,:) = 0;
f_tot(3,:) = 0;
f_th = cross(f_tot,r_all);
f_th = f_th(3,:);
tens = sum(f_th .* th_all) / 2 / pi