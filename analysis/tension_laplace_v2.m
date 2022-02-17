%% Calculate tension by centripetial forces on the anchors
if mod(length(fanc),3) ~=0
    error('the length of fc is not a multiple of 3')
end
fc3 = vec2mat(fc,3)';

rhat = zeros(size(fc3));
r3 = [rbead,rmyo,rmyp];
for i = 1:size(rhat,2)
    rhat(:,i) = r3(:,i) / norm(r3(:,i));
end

centrifugal_f = - sum(fc3 .* rhat);
for i = 1:numel(centrifugal_f)
    if norm(r3(1:2,i)) > r_ring
        centrifugal_f(i) = centrifugal_f(i) ;%- kwall * (norm(r3(1:2,i))-r_ring);
    end
end
tens = sum(centrifugal_f)/2/pi;