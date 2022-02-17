function [theta3,r3,z3,flag]=cutWhiskers(theta4,r4,z4,r_ring)
	whiskerFlagR = (r4 <= (median(r4)-0.2));
	whiskerFlagZ = (z4-median(z4)).^2 >= 0.2^2;
	flag = ~(whiskerFlagR | whiskerFlagZ);
	r3 = r4(flag);
	theta3 = theta4(flag);
	z3 = z4(flag);
end