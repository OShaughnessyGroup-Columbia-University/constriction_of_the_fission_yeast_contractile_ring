tau = 2/30.;
ga_myo2 = 40;
ga_myp2 = 40;

alphafor = 100; % formin kinetic scheme parameter i.e. k1/koff
epsilonfor = 0.0012; % formin kinetic scheme parameter i.e. k2/k1
	% with these parameters, formin:node ratio is 1.1 and
	% ratio of 1 to 2 nodes at steady state is 0.1, ratio
	% of 1 to zero nodes is 0.01.

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
    fhead = 1.75;
    fheadmyp = 1.;
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
    myo2p_unbinding_force = 40;
    kcap = myo2p_unbinding_force / (rcapmyo_short-r_flat_myo);
    r_flat_myp = .5*rcapmyp;
    myp2p_unbinding_force = 30;
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
    rexc = 2*rcapmyp; % careful: this is actually the excluded "diameter"! Not radius!
    rcdc15 = 0.035;
    rexc_n = rcapmyo_long*2;
    rexc_po = rcapmyp + 2*rcapmyo_short;
    kexc_oo = 40*3;
    kexc_pp = 320;
    max_repulsion_force = 6;
    kexc_po = 100;
    % excluded volume for nodes
    % actin excluded volume
    kex_act = 10;
    rex_act = 0.015;
    qmax = pi/6.;
    %% constraints
    % septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
    vs_abs = 0*.07/60;
    % crosslinker spring constant (25 in Bala. for a-actinin)
    kx = 25;
    % corsslinker rest length
    rx0 = 0.03;
    
    % Myosin cluster drag coefficient from membrane (1300 in exp)
    gm = 500; % best value 777 under amplification factor 2
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
        %gb = 8*pi*half_length*mu/(sigma-gamma);%=0.28
        gb = 0.2*(3.1); % adjustment of joseph's value
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
