% change from parameter_std: added myosin atp cycle time
% change from v1: added node excluded volume according to FPALM Cdc15p
% change from v2: no ga
% change from v3: set node excluded volume to be the same as Myo2p size. if
% a smaller excluded volume according to Cdc15p size is used, nodes
% aggregate too much, even if a high membrane drag is used such that Myo2p
% speed is lower than experimentally seen. Also, r_flat_myo and r_flat_myp
% are defined. 
% change from v4: use a separate kexc_po that's scanned
% what happened to v5?
% change from v6: calculate r_flat_myo and r_flat_myp as half of capture radii. 3d scan for rcapmyp, kcap_p and kexc_po
% v8-v10 skipped
% change from v7: calculate kcap and kcap_p by energy
% change from v11: calculate kexc_po by energy, and scale gm and gf
% change from v12: 3d scan binding energy of Myo2p and actin, Myp2p and
% actin, and Myo2p-Myp2p repulsion energy
% change from v13: 2d scan of Myo2p-actin binding energy and actin excluded
% volume
% change from v14: 2d scan of Myp2p-actin binding energy and actin excluded
% volume
% change from v15: no scan, use best values
% change from v16: use force language for Myo2p-actin binding, Myp2p-actin
% binding and Myo2p-Myp2p repulsion. scan myp2p_unbinding_force and max_repulsion_force
% change from v17: no scan, use best values for new rule
% change from v18: fhead, fone, unbinding forces and repulsion forces X2
% v19 to v27 skipped
% change from v28: 2d scan of anchor drags to reproduce experimental node
% velocites 
% v29 and v30 skipped
% change from v31: rexc between nodes made larger
% v32 to v34 skipped
% change from v35: use best value of gm (777). 2d scan of
% myo2p_unbinding_force and myp2p_unbinding_force
% change from v36: use 2.1 pN for fmyo and fmyp, lowered the scan range to
% 40 pN max
% v37 and v38 skipped 
% change from v39: use 20 pN for myo2 and myp2 unbinding forces. Scan actin
% drag and myp2 drag.
% change from v40: no scan, use max drag in v40
% v41 skipped
% change from v42: do not define atp_cycle
% v43 to v50 skipped
% change from v51: larger kwall
% change from v52: change tau_o_dt
% change from v53: max_repulsion_force is set to zero (used to be 6). This
% sets myo2-myp2 repulsion.
% change from v54: unbinding forces = 60 pN
% change from v55: use myp2 unbinding force 10 pN, gm 2000
% v56 skipped
% change from v57: scan myo2 unbinding force
%% scan myo2p_unbinding_force
x1 = 2:2:64; 
x2 = 1:8;
[X1,X2] = meshgrid(x1,x2);

myo2p_unbinding_force = X1(itask);
%% Constants
% Misc.
    % tuning coefficient for constraint feedback
    alpha = -10;         
    if ~exist('nsp', 'var')
        nsp = 1;
    end
    imax = 1e5;          % for debug. can be deleted later.
    tau_o_dt = 2;            % in multiples of dt
% Forces
    % Static force of a myosin head on a filament (pN) (1 in exp)
    fhead = 2.1;
    fheadmyp = 2.1;
    % number of heads per myosin cluster
    nhead = 40 / nsp;
    % maximal stall force per myosin cluster per filament
    fone = 4*2;
    % Maximum number of filaments that a myosin cluster interacts with, without overloading (10 in dev cell)
    nsat = fhead*nhead/fone;
    % load-free velocity
    vmyo0 = 0.24;
    % Km and vm of Michaelis-Menton equation, vmyo0 vs number of heads per
    % filament (11.3 and 0.35 in Stark 2010)
    km = 11.3;
    vm = 0.35;
    % spring constant for myosin - actin capture force
    rcapmyo_long = .132/2;
    rcapmyo_short = .102/2;
    rcapmyo_u = 0.1; % capture radius of (spherical) unanchored myo2 clusters. Not used unless in det ring sim.
    rcapmyp = 0.1;
    r_flat_myo = .5*rcapmyo_short;
    r_flat_myp = .5*rcapmyp;
%     myo2p_unbinding_force = 60;
    myp2p_unbinding_force = 10;
    kcap = myo2p_unbinding_force / (rcapmyo_short-r_flat_myo);
    kcap_p = myp2p_unbinding_force / (rcapmyp-r_flat_myp);
    
    select4 = false;     % if true, myp2 interacts with 4 actin filaments only

    % membrane elastic constant
    kwall = 1000;
    % breaking length of xlinkers (.05 in dev cell)
    l_break = .05;
        l_break_sq = l_break^2;    
    % peeling force
    fpeel = 1e8;
    % myosin exluded volume
    rexc = rcapmyp + rcapmyo_short;           % careful: this is actually the excluded "diameter"! Not radius!
    rexc_n = rcapmyo_long*2;
    kexc = 40*4;
    max_repulsion_force = 0*6;
    kexc_po = max_repulsion_force / rexc;
    % excluded volume for nodes
    % actin excluded volume
    rex_act = 0.02;
    kex_act = 20*2;
    %% constraints
    % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
    vs_abs = 0*.07/60;
    % crosslinker spring constant (25 in Bala. for a-actinin)
    kx = 25;
    % corsslinker rest length
    rx0 = 0.03;
    
    % Myosin cluster drag coefficient from membrane (1300 in exp)
    gm = 2000; % best value 777 under amplification factor 2
    % Formin drag coefficient (1900 in exp)
    gf = 0;
        % water viscosity (in pN um-2 / s)
        mu = 0.6;
        % solution drag of actin beads, see Broersma 1960
        half_length = dbead / 2;
        dact = 0.010;            % diameter of actin filament
        b = dact / 2;
        sigma = log(2*half_length / b);
        gamma = 0.35 - 4 * (1/sigma - 0.43)^2;
        gb = 3.1;
        % solution drag of myosin clusters 
        gmsol = 78;
        % bending stiffness of actin filament
        kappa = 0.041;
        % kappa over dbead^2
        kod2 = kappa / (dbead^2);
        % maximal bending force on a bead
        maxfb = 1e6; % 0.14 at 3kT
% turnover
    % Actin severing rate per filament length by cofilin (46% actin left
    % after 60s of constriction. cf. 1.8 / 60 um-1 s-1 in dev cell)
    rsev = .93/60;
    % off rate of formins from nodes (zero in our latest model. if this is nonzero, 
    % remember to change binding_rate_for)
    kofffor = 0*0.023; 
    % take-out rate of nodes (mean of myosin and formin off rates in dev cell, 0.023 and 0.026 respectively.)
    koffmyo = 0.0245;
    koffmyp = 0.026;
    % take-out rate of a-actinin (3.3 in Dev Cell.)
    koffx = 3.3;
    % binding rates
    binding_rate_myo2 = koffmyo * rhom; %(0.2 in dev cell)
    binding_rate_myo2_init = binding_rate_myo2;
    binding_rate_myp2 = koffmyp * rhom * 2/3;
    binding_rate_myp2_init = binding_rate_myp2;
    binding_rate_for = koffmyo * rhof; % 0.35 in dev cell
    binding_rate_for_init = binding_rate_for;
    binding_rate_x = koffx * rhoain;
    vpol = 0.127; % (0.07 in dev cell)
    vpol_init = vpol;
    
    dbead_first = repmat(dbead,1,numel(ifor));
    
% initial velocity = 0
    vbead = 0 * rbead;
    vmyo = 0 * rmyo;
    fc = zeros(3*(size(rbead,2)+size(rmyo,2)),1);
    poly_time = false;