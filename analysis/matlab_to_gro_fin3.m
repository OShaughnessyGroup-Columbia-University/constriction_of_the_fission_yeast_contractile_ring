% Script to convert a .mat file with variables rbead, ifor, ipt, rmyp and rmyo...
% to a GROMACs file to be used by VMD. Output are 2 text files with gromacs (.gro)
% format and psf format. Please see example gro and psf files to understand formatting 

%filename is the input .mat file as the function argument
function fname_gro = matlab_to_gro_fin3(filename)

% loading the file   
    load(filename)
    xmat(1,2)=1;
    xmat(2,1)=1;
    
    %manipulating the variable rbead to obtain xyz coordinates of the actin filaments
    len = length(rbead(1,:));
    
    %PARAMETERS
    
    lvl1 = 8; %myo2 # of branches in upper (long) branch level 
    lvl2 = 2; % myo2 # of branches in bottom (short) branch level
    
    %index column for actin beads for the gromacs file
    ind = (1:len);
    ind = ind.';
    
    %actin 'bead' xyz coordinate columns
    if exist('rbead')==1 && isempty(rbead) == 0
        rbeadtrans = rbead.';
    else
        rbeadtrans = 0;
    end
    
    rbeadtrans = rbeadtrans.*100; %scaling all coordinates by 100 to make VMD visualization better
    indnew = [ind rbeadtrans];
    
    %%-------------------------------------------------------------------------------

    %GROMACS output file initialization
    x = filename(1:end-4);
    x = strcat(x,'.gro');
    
    fileID = fopen(x, 'w');
    fname_gro = fileID;
    
    %formatting output file according to gromacs specifications
    fmt='%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n';
    
    fprintf(fileID,strcat(filename,'  Ring simulation, t= 0.0\n'));
    
    %length of rmyo
    rmyot = rmyo.';
    myol = length(rmyot(:,1));
    
    %length of crosslinks
    [q, w] = find(tril(xmat));
    cross = [q;w];
    crosslen = length(cross(:,1));
    
    
    %Total number of atoms is given by variable tot_num. Accurate value needed for gro and psf files to run
    % tot_num is found by addding all the various parameter values above and from the .mat variables
    if exist('rmyp') ==1
        tot_num = (lvl1*3*myol) + (2*myol) + (lvl1*2*myol)+ len +...
			length(ifor)+ ((32)*length(rmyp(1,:))) + crosslen;
    else
        tot_num = (lvl1*3*myol) + (2*myol) + (lvl1*2*myol) + len + length(ifor) + crosslen;
    end
    
    fprintf(fileID,'%5d\n',tot_num);
        
    % -------------------------------------------------------------------
    
    % ACTIN section
    %column 1 bead indices
    mol_ind(1:len) = 1;
    
    l = 1; %counting variable
    for j = 1:len
        if j <= length(ifor)
          mol_ind(ifor(j) : ipt(j)) = l;
          l = l+1;
          j = j+1;
        
        end
    end
    

    mol(1:len,2) = 'A';
    %column 3 bead label for the ACTIN beads (formin and pointed end added as ghost beads later)
    for i = 1:len
        mol(i,:) = 'CA';
    end
    
    %GROMACS body printing all columns
    for i= 1:len   
        fprintf(fileID, fmt,mol_ind(i),'ACTIN',strcat('    ',mol(i,:)),indnew(i,:)); 
    end

    %---------------------------------------------------------------------
    
    %MYP2 section:
    
    if exist('rmyp') == 1 && isempty(rmyp) == 0 %to account for runs without myp2 present

        
        %transpose of rmyp 
        rmyptrans = rmyp.';
        rmyptrans = rmyptrans.*100;
        myplen = length(rmyptrans(:,1));
        %making matrix for myp2 points
        myp_points((32*myplen),3) = 0;
        myp_points_len = length(myp_points(:,1));
        
        
        %inserting 8 dimers at correct locations, all 64 points generated by myp2_shape.m
        % added to the variable myp_points. for details of creation see myp2_shape file 
		for i = 1:2:size(rmyptrans,1)
            		myp_points(((i-1)*32+1):(((i-1)*32+1)+63),:) = myp2_shape(rmyptrans(i,1),...
                	rmyptrans(i,2),rmyptrans(i,3),rmyptrans(i+1,1)...
                	,rmyptrans(i+1,2),rmyptrans(i+1,3));
        	end
      
        %GROMACS file info for myp2
        
        %indexing
        m_ind = ((length(ind)+1):((length(ind)+myp_points_len)));
        m_ind = m_ind.';
        
        %identity of each bead to corresponding dimer (first 64 beads to first dimer and so on)
        m_id = ((len+1):myp_points_len);
        m_id = m_id.';

       for i = 1:64:myp_points_len
           m_id(i:(i+63))= ceil(i/64);       
       end

		
        
       %combining indices and coordinates
        m_ind_r= [m_ind myp_points];

        %naming: myp2 beads are given 4 names. details of names can be found in the ActoMyosin.tcl file 
        
        %column for storing names 
        myp_name(1:myp_points_len,2) = 'A';
        
        for i = 1:myp_points_len
            myp_name(i,:)= 'CM';
        end

        for i = 1:myp_points_len
            if mod(i,8)== 1 || mod(i,8)== 2  
                myp_name(i,:) = 'CM'; 
            
            elseif mod(i,8)== 3 || mod(i,8)== 6
                myp_name(i,:) = 'CM';
                
            elseif mod(i,8)== 4 || mod(i,8)== 7
                myp_name(i,:) = 'CE';
                
            elseif mod(i,8)== 5 || mod(i,8)== 0
                myp_name(i,:) = 'CH';
            end
            
        end
        
       %Converting the myp_points coordinates into gromacs specified columns 
       myp_bead(myp_points_len,3) = 0;
       
       for i = 1:myp_points_len
           if mod(i,8) == 3 || mod(i,8) == 6
               myp_bead(i,:) = myp_points(i,:);
           end
       end
       
       mb = myp_bead(:,1) == 0;
       myp_bead(mb,:) = [];
       

        
    %gromacs print statement
        for i = 1:size(m_ind_r(:,1))
            fprintf(fileID, fmt,m_id(i),'MYPII',strcat('    ',myp_name(i,:)),m_ind_r(i,:));   
        end
        
    else
        myp_points_len = 0;
        myp_points = 0;
        m_ind_r = 0;
    end
        
   %---------------------------------------------------------------------    

    %MYO2  BOUQUET STRUCTURES
    
    lvl1 = 8; %same as above in parameter section
    lvl2 = 2; %same as above
    
    lvl1_len = 0.096*100; %long branch length is 96 nm, scaled appropriately
    lvl2_len = 0.02*100;  %short branch length is 20nm, scaled appropriately 
    
    % generate p_sept if absent, p_sept is needed for positioning the bouquet against septum
    d_for = 0.044;
    d_myo = 0.094;
    if ~exist('p_sept','var')
        p_sept=[rbead(:,ifor)*r_ring/(r_ring-d_for),rmyo*r_ring/(r_ring-d_myo)];
    end

    %coordinates of the bouquet anchor 'heads'
    p_sept_trans = p_sept.';
    p_sept_myo = p_sept_trans(((length(ifor)+1):end),:);
    rmyotrans = p_sept_myo;
    rmyotrans = rmyotrans.*100; %scaling
    
    myolen = length(rmyotrans(:,1));
    
    rmyo_og = rmyo.';
    
    %matrix containing coordinates of rmyo structures generating using gen_bouquet
    bouquet(((lvl1*3*myolen)+myolen),3) = 0;
    
    %populating the bouquet matrix, first step is adding the 'anchor' bead
    bp = lvl1*3;
    for i = 1:length(bouquet(:,1))
        if mod(i,(bp+1)) == 1
            bouquet(i,:) = rmyotrans(ceil(i/(bp+1)),:); 
        end        
    end

    %adding rest of the points using gen_bouquet
    for i = 2:(bp+1):length(bouquet(:,1))
        bouquet(i:(i+bp-1),:)= gen_bouquet(rmyotrans(ceil(i/(bp+1)),1),...
            rmyotrans(ceil(i/(bp+1)),2), rmyotrans(ceil(i/(bp+1)),3),...
            lvl1_len, lvl2_len, lvl1, lvl2,rmyo_og(ceil(i/(bp+1)),1),...
            rmyo_og(ceil(i/(bp+1)),1),rmyo_og(ceil(i/(bp+1)),1));
    end

    %making the remaing gromacs columns for bead indices and names
    myo_ind = ((len+(myp_points_len)+1):...
        (len+myp_points_len+length(bouquet(:,1))));
    myo_ind = myo_ind.';

    
   %column to ascribe each bead to respective bouquet superstructure
    myo_id = (1:length(bouquet(:,1)));
    myo_id = myo_id.';
    
    for i = 1:(bp+1):(length(bouquet(:,1))-1)
        myo_id(i:(i+bp+1)) = ceil(i/(bp+1));
    end
    
    %column for name of each bead, all head beads given name N (these are ghost beads to preserve bond thickness)
    %except anchor beads, given R and middle branch given Q
    myo_name(1:length(bouquet(:,1)),2)= 'N';
    for i = 1:length(bouquet(:,1))
        myo_name(i,:) = 'CN';
    end
    
    for i = 1:length(myo_ind(:,1))
        if mod(i,(bp+1)) == 1
            myo_name(i,:) = 'CR';
        end       
    end
  
    for i = 2:(bp+1):size(myo_ind(:,1))
        myo_name((i:(i+lvl1-1)),:) = ['CQ';'CQ';'CQ';'CQ';'CQ';'CQ';'CQ';'CQ'];
    end
    
    %-----------------------------------------------------------------
    %adding separate actual beads for heads with name CN
    myo_head(length(bouquet(:,1)),3) = 0;
    

    for i = 1:length(bouquet(:,1))
        if myo_name(i,:) == 'CN'
            myo_head(i,:) = bouquet(i,:);
        end
    end
 
   
    mt = myo_head(:,1) == 0;
    myo_head(mt,:) = [];
   
    
    %=------------------------------------------------------------------
    % remove all names because all are ghost beads for now (once positions are set). Will add separate ghost beads later
    
    for i = 1:length(bouquet(:,1))
        myo_name(i,:) = 'CQ';
    end
    
    %--------------------------------------------------------------------
    
    
    
   %gromacs column printing
    myo_ind_r = [myo_ind bouquet];
   
    for i = 1:size(myo_ind_r(:,1))
       fprintf(fileID, fmt,myo_id(i),'MYOII',strcat('    ',myo_name(i,:)),myo_ind_r(i,:));   
    end
    
    %----------------------------------------------------------------------
    %SHADOW BEADS: used because VMD allows same color and thickness of bond for only 
    %homogenous bonds and hence these ghost or shadow beads needed to control radius and color

    %replacement for shadow beads for MYO2 anchors
    shadow_1 = myo_ind(end);
    %indexing
    anchor = (shadow_1+1 : shadow_1+myolen);
    anchor = anchor.';
    anc_ind = (1:myolen);
    anc_ind = anc_ind.';
    %assignment of each bead and name
    anc_id(1:myolen,2) = 'A';
    for i = 1:myolen
        anc_id(i,:) = 'CB';
    end
    
    %gromacs column printing with coordinates
    p_sept_myo = p_sept_myo*100;
    myo_anc = [anchor p_sept_myo];
    
    for i = 1:size(myo_anc(:,1))
       fprintf(fileID, fmt,anc_ind(i),'MYOII',strcat('    ',anc_id(i,:)),myo_anc(i,:));   
    end
    
    
    %crosslink beads between actin filaments
    
    %xmat is sparsely populated matrix giving bond location of crosslinking between actin
    [q, w] = find(tril(xmat)); %finding all non-empty bonds
    cross = [q;w];
    %populating crosslink matrix with coordinates
    if ~isempty(cross)
        crosslink(length(cross(:,1)),3) = 0;
        for i = 1:length(cross(:,1))
            crosslink(i,:) = rbeadtrans(cross(i),:);
        end
    else
        crosslink = [];
    end
    %indexing, assignments and names as per gromacs requirements    
    link = anchor(end);
    cross_i = 1:crosslen;
    cross_i = cross_i.';
    cross_ind = (link+1 : link+crosslen);
    cross_ind = cross_ind.';
    cross_main = [cross_ind crosslink];
    
    %gromacs print statement
    for i = 1:crosslen
        fprintf(fileID, fmt,cross_i(i),'CROSS','   CX',cross_main(i,:));
    end
    
    %crosslinker bond information to be used in psf file later
    cross_ind1 = cross_ind(1):cross_ind(crosslen/2);
    cross_ind1 = cross_ind1.';
    cross_ind2 = cross_ind((crosslen/2)+1):cross_ind(end);
    cross_ind2 = cross_ind2.';
    %matrix containing crosslinked bonds
    crossbond = [ cross_ind1 cross_ind2];
    
    
    %Formin actual beads to control radius and color, named CF
    formin_p(length(ifor),3) = 0;

    %coordinates assignment
    for i = 1:length(formin_p(:,1))
        formin_p(i,:) = rbeadtrans(ifor(i),:);
    end

    %gromacs file indexing, assignment, naming for formin
    form = cross_ind(end);
    form_i = 1:length(ifor);
    form_i = form_i.';
    form_ind = (form+1 : form+length(ifor));
    form_ind = form_ind.';
    form_main = [form_ind formin_p];

    %gromacs print statement
    for i = 1:length(ifor)
        fprintf(fileID, fmt,form_i(i),'FORMI','   CF',form_main(i,:));
    end
    
    %myo2 bouquet heads actual bonds to control color and radius, named CN
    %indexing, assignment and naming
    myoh = form_ind(end);
    myoh_i= 1:length(myo_head(:,1));
    myoh_i = myoh_i.';
    myoh_ind = (myoh+1: myoh+length(myo_head(:,1)));
    myoh_ind = myoh_ind.';
    myoh_main = [myoh_ind myo_head];
    
    %gromacs print statement
    for i = 1:length(myo_head(:,1))
        fprintf(fileID, fmt,myoh_i(i),'MYOHE','   CN',myoh_main(i,:));
    end

    
    %last line for GROMACS to end file
    fprintf(fileID,'   100.0 100.0 100.0\n');
    fclose(fileID);
    
    %---------------------------------------------------------------------
    
    % PSF FILE STARTS HERE

    y = filename(1:end-4);
    y = strcat(y,'.psf');
    
    %open file
    fileIDp = fopen(y, 'w');
    fname_psf = fileIDp;
    
    %psf file format requirements
    fprintf(fileIDp,'PSF CMAP\n');
    fprintf(fileIDp,'\n');
    fprintf(fileIDp,'       6 !NTITLE\n');
    fprintf(fileIDp,' REMARKS original generated structure x-plor psf file\n');
    fprintf(fileIDp,' REMARKS 2 patches were applied to the molecule.\n');
    fprintf(fileIDp,'\n');
    fprintf(fileIDp,'%5d !NATOM\n', tot_num); 
    
    %format
    fmtp = '   %5d %s%3d    %s    %s       %d   %8.6f       %7.4f           %d\n';
    
    
    %Printing (listing) atoms for the PSF file
    
    %actin
    for i= 1:len   
        fprintf(fileIDp, fmtp,ind(i),'ACT',1,mol(i,:),mol(i,:),1,0.000000,1.0000,0); 
    end
    
    %myp2
    if exist('rmyp') == 1 && isempty(rmyp) == 0
        for i = 1:size(m_ind_r(:,1))
            fprintf(fileIDp, fmtp,m_ind(i),'MYP',1,myp_name(i,:),myp_name(i,:),1,0.000000,1.0000,0);      
        end
    end
    
    %myo2
    for i = 1:size(myo_ind_r(:,1))
        fprintf(fileIDp, fmtp,myo_ind(i),'MYO',1,myo_name(i,:),myo_name(i,:),1,0.000000,1.0000,0); 
    end
    
    %myo2 anchor bead
    for i = 1:size(myo_anc(:,1))
        fprintf(fileIDp, fmtp,anchor(i),'MYO',1,'CB','CB',1,0.000000,1.0000,0); 
    end
    
    %crosslinkers
    for i= 1:crosslen
        fprintf(fileIDp, fmtp,cross_ind(i),'CRO',1,'CX','CX',1,0.000000,1.0000,0);
    end
    
    %formin
    for i= 1:length(ifor)
        fprintf(fileIDp, fmtp,form_ind(i),'FOR',1,'CF','CF',1,0.000000,1.0000,0);
    end
    
    %myo2 heads
    for i= 1:length(myo_head(:,1))
        fprintf(fileIDp, fmtp,myoh_ind(i),'MYH',1,'CN','CN',1,0.000000,1.0000,0);
    end
    
    fprintf(fileIDp,'\n');

    %--------------------------------------------------------------------------------------
    %PSF file requires bond information, 4 pairs of atoms giving 4 bonds giving 8 columns per row
    %Hence each 'bond' is a pair of atoms and this information is determined below

    %actin bonds
    moli(1:len,2) = 'A';
    for i = 1:len
        moli(i,:) = 'CA';
    end    
   
    for i= 1:len
        
        %formin
        if ind(i==ifor)
            moli(i,:) = 'CF';
        end
        
        %pointed end
        if ind(i==ipt)
            moli(i,:) = 'CP';
        end                
    end
    
    %Only pointed ends are not bonded to the next atom
    act_bonds((len+ length(crosslink(:,1))),2) = 0;
    for i = 1:len
        if moli(i,2) ~= 'P'
            act_bonds(i,:) = [ind(i), ind(i+1)]; 
        end        
    end
    
    pt = act_bonds(:,1) == 0;

    %list of actin bonds in the act_bonds matrix
    act_bonds(pt,:) = []; %final
    
    %--------------------------------------------------------------------
    
    % myp2 bonds
    
    if exist('rmyp') == 1 && isempty(rmyp)== 0
        mypbond(myp_points_len,2) = 0;
        
    %arranging bonds for myp2 dumbell structures   
    for i = 1:myp_points_len
        if (mod(i,8))==1
            mypbond(i,:) = [m_ind(i), m_ind(i+1)];
            mypbond(i+1,:) = [m_ind(i+1), m_ind(i+2)];
            mypbond(i+2,:) = [m_ind(i+1), m_ind(i+5)];
        end  
                   
   end
    %mypbond contains all the myp2 bonds                
    pt = mypbond(:,1) == 0;
    mypbond(pt,:) = [];       
        
    else
        mypbond = [];
    end
    
    %--------------------------------------------------------------------
    
    %Myo2 BOUQUET BONDS
    bouqbond(lvl1*3*length(rmyo(1,:)),2) = 1;
    
    %LEVEL 1
    p = (1:(bp):size(myo_ind));
    p = p.';
    q = (1:(bp+1):size(myo_ind));
    q = q.';
    r = (2:(bp+1):size(myo_ind));
    r = r.';

    %bouqbond contains all the myo2 bonds
    for i = 1:size(q(:,1))
        bouqbond((p(i):(p(i)+lvl1-1)),1)= myo_ind(q(i));
        bouqbond((p(i):(p(i)+lvl1-1)),2)= myo_ind(r(i)):myo_ind(r(i)+lvl1-1);
    end
    
    
    %LEVEL 2
    % There are a total 16 heads in each bouquet structure, the following assigns bonds to connect branches
    l = (lvl1+1):(lvl1*3):size(myo_ind);
    l = l.';

    
    m = (lvl1+1):2:size(myo_ind);
    m = m.';
    for i = 1:length(r(:,1))
        bouqbond(l(i):(l(i)+1),1) = myo_ind((q(ceil(l(i)/bp)))+1);
        bouqbond((l(i)+2):(l(i)+3),1)= myo_ind((q(ceil(l(i)/bp)))+2);
        bouqbond((l(i)+4):(l(i)+5),1) = myo_ind((q(ceil(l(i)/bp)))+3);
        bouqbond((l(i)+6):(l(i)+7),1) = myo_ind((q(ceil(l(i)/bp)))+4);
        bouqbond((l(i)+8):(l(i)+9),1) = myo_ind((q(ceil(l(i)/bp)))+5);
        bouqbond((l(i)+10):(l(i)+11),1) = myo_ind((q(ceil(l(i)/bp)))+6);
        bouqbond((l(i)+12):(l(i)+13),1) = myo_ind((q(ceil(l(i)/bp)))+7);
        bouqbond((l(i)+14):(l(i)+15),1) = myo_ind((q(ceil(l(i)/bp)))+8);
        
        bouqbond(l(i):(l(i)+1),2) = [myo_ind((q(ceil(l(i)/bp)))+9),myo_ind((q(ceil(l(i)/bp)))+10)];
        bouqbond((l(i)+2):(l(i)+3),2)= [myo_ind((q(ceil(l(i)/bp)))+11),myo_ind((q(ceil(l(i)/bp)))+12)];
        bouqbond((l(i)+4):(l(i)+5),2) = [myo_ind((q(ceil(l(i)/bp)))+13),myo_ind((q(ceil(l(i)/bp)))+14)];
        bouqbond((l(i)+6):(l(i)+7),2) = [myo_ind((q(ceil(l(i)/bp)))+15),myo_ind((q(ceil(l(i)/bp)))+16)];
        bouqbond((l(i)+8):(l(i)+9),2) = [myo_ind((q(ceil(l(i)/bp)))+17),myo_ind((q(ceil(l(i)/bp)))+18)];
        bouqbond((l(i)+10):(l(i)+11),2) = [myo_ind((q(ceil(l(i)/bp)))+19),myo_ind((q(ceil(l(i)/bp)))+20)];
        bouqbond((l(i)+12):(l(i)+13),2) = [myo_ind((q(ceil(l(i)/bp)))+21),myo_ind((q(ceil(l(i)/bp)))+22)];
        bouqbond((l(i)+14):(l(i)+15),2) = [myo_ind((q(ceil(l(i)/bp)))+23),myo_ind((q(ceil(l(i)/bp)))+24)];
    end
   
    % -------------------------------------------------------------------
    
    %act_bonds IS THE MAIN MATRIX THAT COLLECTS ALL THE BOND INFO TOGETHER:
       
    act_bonds= [act_bonds; mypbond; bouqbond; crossbond];
    
    % the following procedure cuts the matrix repeatedly to make 8 columns in each row with 1's at the very end to fill matrix
    if mod(size(act_bonds(:,1)),2) == 1
        cut = (size(act_bonds(:,1))+1);
        ab1 = act_bonds(1:(cut/2), :);
        ab2 = act_bonds(((cut/2)+1):end, :);
        f1 = [1, 1];
        ab2 = [ab2 ; f1];
        ab = [ab1 ab2];       
    else
        cut = size(act_bonds(:,1));
        ab1 = act_bonds(1:(cut/2), :);
        ab2 = act_bonds(((cut/2)+1):end, :);
        ab = [ab1 ab2]; 
    end
    
    if mod(size(ab(:,1)),2) == 1
        cut2 = (size(ab(:,1))+1)/2;
        ab12= ab(1:cut2, :);
        ab22 = ab((cut2+1):end,:);
        fill = [1,1,1,1];
        ab22 = [ab22 ; fill];
        ab2 = [ab12 ab22];
    else
        cut2 = (size(ab(:,1)))/2;
        ab12= ab(1:cut2, :);
        ab22 = ab((cut2+1):end,:);
        ab2 = [ab12 ab22];
    end
    
    fmtb= '    %4d    %4d      %4d      %4d      %4d      %4d      %4d      %4d\n';
    
    actbonds = (size(act_bonds(:,1)));
    
    
    % PRINT statement FOR THE PSF FILE
    fprintf(fileIDp, '   %5d !NBOND: bonds\n', actbonds(1));
    for i = 1:size(ab2(:,1))
        fprintf(fileIDp, fmtb, ab2(i,:));
    end
    %---------------------------------------------------------------------------------------------------------------------
    %The septum.m file takes data from the .mat file to generate a .tcl file with the triangular coordinates for the septum  
    septum(filename)
   
   
end


