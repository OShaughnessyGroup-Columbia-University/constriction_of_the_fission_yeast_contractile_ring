function good = isgood(ifil, rbead, rmyo, ifor, ipt, r_ring)
good = 1;
curr_flt = rbead(:,ifor(ifil):ipt(ifil));

% extract velocities and myosin positions;
% rbead(3,:) = rbead(3,:) - mean(rbead(3,:));
% rmyo(3,:) = rmyo(3,:) - mean(rmyo(3,:));
% 
% [theta,rho,z] = cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
% theta = theta + pi/2;
% dirn = [cos(theta);sin(theta);0*(theta)];
% % if filament is a whisker
% if(any(abs(curr_flt(3,:)) > 100/1000))
%     good = 0;
% end
% if(any(sqrt(curr_flt(1,:).^2 + curr_flt(2,:).^2) < r_ring - 200/1000))
%     good = 0;
% end

% if filament is too short
% if(size(curr_flt,2) < 7) % filament must have at least n beads
%     good = 0;
% end

end