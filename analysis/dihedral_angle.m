function psi = dihedral_angle(a, b, c, d)
ba = b - a;
cb = c - b;
dc = d - c;

n1 = cross(ba, cb);
n1 = n1/norm(n1);
n2 = cross(cb, dc);
n2 = n2/norm(n2);
m1 = cross(n1, cb)/norm(cb);

x = sum(n1.*n2);
y = sum(m1.*n2);

psi = atan2(y, x);