inside_outside_ellipse_myo;
reachFlagMyo = reachFlag;
naiveProjectedReachFlagMyo = naiveProjectedReachFlag;
figure;
inside_outside_ellipse_myp;
reachFlagMyp = reachFlag;
naiveProjectedReachFlagMyp = naiveProjectedReachFlag;

nBeads = size(r3);

atleast_1_myo_1_myp = size(r3(1,reachFlagMyo > 0 | reachFlagMyp > 0 ))/nBeads * 100
atleast_1_myo_1_myp_intg_theta = size(r3(1,naiveProjectedReachFlagMyo > 0 | naiveProjectedReachFlagMyp > 0 ))/nBeads * 100

no_myo_myp = size(r3(1,reachFlagMyo == 0 & reachFlagMyp == 0 ))/nBeads * 100
no_myo_myp_intg_theta = size(r3(1,naiveProjectedReachFlagMyo == 0 & naiveProjectedReachFlagMyp == 0 ))/nBeads * 100

only_myo = size(r3(1,reachFlagMyo > 0 & reachFlagMyp == 0 ))/nBeads * 100
only_myp = size(r3(1,reachFlagMyo == 0 & reachFlagMyp > 0 ))/nBeads * 100
only_myo_intg_theta = size(r3(1,naiveProjectedReachFlagMyo > 0 & naiveProjectedReachFlagMyp == 0 ))/nBeads * 100
only_myp_intg_theta = size(r3(1,naiveProjectedReachFlagMyo == 0 & naiveProjectedReachFlagMyp > 0 ))/nBeads * 100

both_myo_myp = size(r3(1,reachFlagMyo > 0 & reachFlagMyp > 0 ))/nBeads * 100
atleast_1_myo = size(r3(1,reachFlagMyo > 0))/nBeads * 100
atleast_1_myp = size(r3(1,reachFlagMyp > 0))/nBeads * 100

both_myo_myp_intg_theta = size(r3(1,naiveProjectedReachFlagMyo > 0 & naiveProjectedReachFlagMyp > 0 ))/nBeads * 100
atleast_1_myo_intg_theta = size(r3(1,naiveProjectedReachFlagMyo > 0))/nBeads * 100
atleast_1_myp_intg_theta = size(r3(1,naiveProjectedReachFlagMyp > 0))/nBeads * 100