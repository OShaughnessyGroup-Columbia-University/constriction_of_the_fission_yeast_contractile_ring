% change from original: use r = 3.7/2
% change from v1: rhof = 17.2, rhom = 7.5 * 2.5
% v2 to v8 skipped
%% Constants

% Geometry
    % radius of the intact cell (3.7 / 2 in septation paper)
    r = 3.7/2;    
    
    z_up = 0; % distance from the plane of the ring to the equator
%         new_r = sqrt(r^2 - z_up^2);
%         lr = new_r * 2 * pi; % length of the ring
    lr = r * 2 * pi;    
    lf = 1.3; % mean length of actin filament in the ring
    % Initial anchored region(s)
    % la is a 2 by n matrix, each column specifies the start and end of an
    % anchored region
    la = lr * [0; 1];
    
    % binding radius of ain1
    rxbind = 0.05;
    % width of the myo2 ring
    wr = 0.1;
    % depths of formin and Myo2 clusters under the membrane
    d_for = .044;        % 0.044 micron from sup res paper Fig. 2F
    d_myo = .094;        % 0.094 micron from sup res paper Fig. 2D
    d_stem = 0.13;
    d_cdc15 = 0.040;
    % distance between beads on an actin filament
    dbead = 0.1;
        dbead_sq = dbead ^2;
    % capture radius of myosin clusters (unfortunately the name rmyo is
    % used by another variable)
    % minimum distance between myosin clusters (excluded volume)
    rexc = 0;
    
% density of components
    nsp = 2.5;
    rhom = 7.5 * nsp;
    rhof = 17.2;
    rhoain = 25;
    
for itrial = 1:1e5
    
    %% Initial positions of all molecules
    [rbead, rmyo, ipt, ifor, bancf, bancm, xmat, fmmat] = initial_circle_rand (rhom, rhof, rhoain, r, la, wr, d_for, d_myo, dbead, rexc, rxbind);
        rmyp = double.empty(3,0);
        
%         rbead = double.empty(3,0);
%         rmyo = double.empty(3,0);
%         ifor = int16.empty(1,0);
%         ipt = int16.empty(1,0);
%         bancf = logical.empty(1,0);
%         bancm = logical.empty(1,0);
%         xmat = logical.empty(0,0);
%         fmmat = logical.empty(0,0);
        fil_tag = 1:numel(ipt);
        r_ring = r; 
    %% test if the entire ring is bundled
        confbad = false;
        for i = 1:size(rbead,2)
            % find out possible neighbors whose x and y positions are within
            % half of dbead
            xtemp = abs(rbead(1,i)-rbead(1,:)) < .5 * dbead;
            ytemp = abs(rbead(2,i)-rbead(2,:)) < .5 * dbead;
            xytemp = and(xtemp, ytemp);
            xytemp(i) = false;
            % if there is a bead without neighbors, this configuration is bad
            if ~any(xytemp)
                confbad = false%true;
                break
            end
        end
        % test if any two myosin clusters are within excluded radius
        for i = 1:size(rmyo,2)
            xtemp = abs(rmyo(1,i)-rmyo(1,:)) < rexc;    % myosins that are within rexc of the current myosin in the x direction
            ytemp = abs(rmyo(2,i)-rmyo(2,:)) < rexc;
            ztemp = abs(rmyo(3,i)-rmyo(3,:)) < rexc;
            xytemp = and(xtemp, ytemp);
            xyztemp = and(xytemp,ztemp);
            xyztemp(i) = false;                         % now xyztemp are myosins within rexc of the current myosin in x,y and z directions
            if any(xyztemp)     % if there is any such candidate
                for j = find(xyztemp)
                    rmyo_diff = rmyo(:,j) - rmyo(:,i);
                    rmyo_diff_sq = rmyo_diff' * rmyo_diff;
                    if rmyo_diff_sq < rexc * rexc
                        confbad = true;
                        break
                    end
                end
            end
            if confbad
                break
            end
        end
        if ~confbad
            break
        end
end

%% if a good configuration is found, output
if ~confbad
    vbead = zeros(size(rbead));
    vmyo = zeros(size(rmyo));
%     ringplot_pick
%     sphereplot(r)
%     view(0,90)
else
%     ringplot_pick
%     sphereplot(r)
%     view(0,90)
    error('all trials gave bad configurations')
end

nm = size(rmyo,2);                              % number of myosin clusters
rmyp = double.empty(3,0);
vmyp = double.empty(3,0);