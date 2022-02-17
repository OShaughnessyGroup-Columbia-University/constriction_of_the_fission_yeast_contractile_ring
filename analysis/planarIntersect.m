function cross = planarIntersect2D(p0, p1, q0, q1)
% Detect if line segments going from q0 to q1 and p0 to p1 intersect. For
% there to be an intersection, initial and final points of both segments
% must be on opposite sides of the line defined by the other segment.

%% Does q cross the line defined by p
pslope = (p1(2) - p0(2))/(p1(1) - p0(1));
qcrossp = ...
    (q0(2) - p0(2) > pslope*(q0(1) - p0(1)) && q1(2) - p1(2) < pslope*(q1(1) - p1(1))) || ...
    (q0(2) - p0(2) < pslope*(q0(1) - p0(1)) && q1(2) - p1(2) > pslope*(q1(1) - p1(1)));

%% Does p cross the line defined by q
qslope = (q1(2) - q0(2))/(q1(1) - q0(1));
pcrossq = ...
    (p0(2) - q0(2) > qslope*(p0(1) - q0(1)) && p1(2) - q1(2) < qslope*(p1(1) - q1(1))) || ...
    (p0(2) - q0(2) < qslope*(p0(1) - q0(1)) && p1(2) - q1(2) > qslope*(p1(1) - q1(1)));

cross = pcrossq && qcrossp;
end