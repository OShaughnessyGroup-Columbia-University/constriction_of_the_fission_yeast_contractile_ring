function X = planProjCross(p0, p1, q0, q1)
q10 = q1 - q0;
p10 = p1 - p0;
pq0 = p0 - q0;

n = cross(p10, q10);
n = n/norm(n);

%% projection out components along n
p0 = p0 - n*sum(n.*p0);
p1 = p1 - n*sum(n.*p1);
q0 = q0 - n*sum(n.*q0);
q1 = q1 - n*sum(n.*q1);

%% switching to phat, n2 basis
phat = (p1 - p0)/norm(p1 - p0);
n2 = cross(n, phat);

p0 = [sum(p0.*phat), sum(p0.*n2)];
p1 = [sum(p1.*phat), sum(p1.*n2)];
q0 = [sum(q0.*phat), sum(q0.*n2)];
q1 = [sum(q1.*phat), sum(q1.*n2)];

X = planarIntersect2D(p0, p1, q0, q1);

% if X==1
%     figure
%     hold on
%     p = [p0; p1];
%     plot(p(:, 1), p(:, 2), 'ro-' )
%     q = [q0; q1];
%     plot(q(:, 1), q(:, 2), 'ko-' )
% end
end