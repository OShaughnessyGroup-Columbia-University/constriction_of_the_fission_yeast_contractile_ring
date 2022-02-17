function surfPlots3(prf)
    load([prf '_1_00min.mat'])
    mat_files = dir([sprintf('%s_*_%dmin.mat', prf, 10)]);
    nmax = length(mat_files);
    % mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','myo2e2', 10)]);
    % mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','delmyp7', 10)]);
    % mat_files = dir([workingDir '/' sprintf('%s_*_%dmin.mat','sizep5', 10)]);

    clear myo_t myp_t act_t rwhisk 
    % myo_w = zeros(1,nmax);
    % myo_t = zeros(1,nmax);
    % myp_w = zeros(1,nmax);
    % myp_t = zeros(1,nmax);
    % act_w = zeros(1,nmax);
    % act_t = zeros(1,nmax);
    % rwhisk = zeros(1,nmax);

    for i=1:nmax
        load(mat_files(i).name);

        rbead(:,ifor) = []; %#ok<*SAGROW>
        [theta1,r1,z1]=cart2pol(rmyp(1,:),rmyp(2,:),rmyp(3,:));
        [theta2,r2,z2]=cart2pol(rmyo(1,:),rmyo(2,:),rmyo(3,:));
        [theta3,r3,z3]=cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));
%         figure
%         histogram(r2)
%         axis([0.08 0.12 0 length(r2)])

        %% identify whiskers
        dmean = median(r_ring - r3);
        whiskers = (r_ring - r3 > dmean + 2*rcapmyp);
        rwhisk(i) = length(find(whiskers))/length(r3);
        r3(whiskers) = [];
        z3(whiskers) = [];
        whiskers = (r_ring - r1 > dmean + 2*rcapmyp);
%         figure
%         hold on
%         histogram(r1)
%         histogram(r1(whiskers))
        r1(whiskers) = [];
        z1(whiskers) = [];
        rwhisk_myp(i) = length(find(whiskers))/length(r1);
%         histogram(r1)

        %% convolve with tophat
        z = [r1,r2,r3];
        M = zeros(1000*length(z),2);
        t = [ones(1,length(r1)),2*ones(1,length(r2)),3*ones(1,length(r3))];
        wdth = [rcapmyp*ones(1,length(z1)),rcapmyo_short*ones(1,length(z2)),0.0025*ones(1,length(z3))];

        for k=1:length(z)
            M(1000*(k-1)+1:1000*(k-1)+1000,1) = random('Uniform',r_ring-z(k)-wdth(k),r_ring-z(k)+wdth(k),[1000,1]);
            M(1000*(k-1)+1:1000*(k-1)+1000,2) = t(k);
        end

        %% thickness
        myp_t(i) = quantile(M(M(:,2)==1), 0.95) - quantile(M(M(:,2)==1), 0.05);
        myo_t(i) = quantile(M(M(:,2)==2), 0.95) - quantile(M(M(:,2)==2), 0.05);
%         if(myo_t(i) > 0.1)
%             i
%             myo_t(i)
%             
%             figure
%             histogram(M(M(:,2)==2))
%         end
        act_t(i) = quantile(M(M(:,2)==3), 0.95) - quantile(M(M(:,2)==3), 0.05);
        
        %% width
        myp_w(i) = quantile(z1, 0.95) - quantile(z1, 0.05);
        myo_w(i) = quantile(z2, 0.95) - quantile(z2, 0.05);
        act_w(i) = quantile(z3, 0.95) - quantile(z3, 0.05);
    end
    myptave = mean(myp_t)
    myotave = mean(myo_t)
    acttave = mean(act_t)
    
    rwhiskave = mean(rwhisk)
    rwhisk_myp_ave = mean(rwhisk_myp)
    
    mypwave = mean(myp_w)
    myowave = mean(myo_w)
    actwave = mean(act_w)
end