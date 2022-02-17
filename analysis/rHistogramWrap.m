bndl_act_frac = [];
bndl_act_myo_frac = [];
bndl_act_myp_frac = [];
bndl_myp_frac = [];
thick = [];
width = [];
rbndl = [];
nx_list = [];
xsect_fil = [];
myo_myp_sep = [];

dir = '.';
% tag = 'wtfc4';
% nmax = 200;
% tag = 'weak_myp_unbind';
% tag = 'allmyo3';
% tag = 'noga';
% tag = 'rxkx5';
% snap_time = 180;
% unit = 'sec';
% tag = 'delmyp7';
% tag = 'myo2e2';
% tag = 'wtfc4';
% snap_time = 1;

% tag = 'bind_2da';
% tag = 'rass'
% tag = 'assmyp'
% tag = 'allmyo3'
% tag = 'matur'
% tag = 'vcont2';
% snap_time = 2400;
% unit = 'sec';

% for k=i0:i0+nval-1
for k=i0:nval:nmax
% for k=6:6:60
% 	if exist(sprintf('%s/%s_%d_10min.mat',tag,tag,k)) == 2
% 		load(sprintf('%s/%s_%d_10min.mat',tag,tag,k))
	if exist(sprintf(['%s/%s_%d_' num2str(snap_time) unit '.mat'],dir,tag,k)) == 2
        try
            load(sprintf(['%s/%s_%d_' num2str(snap_time) unit '.mat'],dir,tag,k))
        catch
            disp('Could not load file, though it exists!')
        end
		rHist_noplot
		%close all;
		%break;
		bndl_act_frac = [bndl_act_frac cnt/size(rbead,2)*100];
		bndl_act_myo_frac = [bndl_act_myo_frac cntMyo/cnt*100];
		bndl_act_myp_frac = [bndl_act_myp_frac cntMyp/cnt*100];
		bndl_myp_frac = [bndl_myp_frac cntMypInBundle/size(rmyp,2)*100];

		% 95 thickness
		nx_list = [nx_list nxfil];
		thick = [thick prctile(rk,99) - prctile(rk,1)];
		width = [width prctile(zk,99) - prctile(zk,1)];
		rbndl = [rbndl rb];
		
		tot_len = 0;
		
		for k = 1:size(rbead,2)-1
			if any(ifor(:) == k + 1)
				continue;
			end
			tot_len = tot_len + sqrt((rbead(1,k)-rbead(1,k+1))^2 + (rbead(2,k)-rbead(2,k+1))^2 + (rbead(3,k)-rbead(3,k+1))^2);
		end
		xsect_fil = [xsect_fil tot_len/(2*pi*r_ring)];
		myo_myp_sep = [myo_myp_sep sep];
    else
        disp(['Could not load file ' tag '_' num2str(k) '_' num2str(snap_time) unit '.mat'])
        nx_list = [nx_list nan];
        thick = [thick nan];
        width = [width nan];
        rbndl = [rbndl nan];
        bndl_act_frac = [bndl_act_frac nan];
        bndl_act_myo_frac = [bndl_act_myo_frac nan];
        bndl_act_myp_frac = [bndl_act_myp_frac nan];
        bndl_myp_frac = [bndl_myp_frac nan];
        xsect_fil = [xsect_fil nan];
        myo_myp_sep = [myo_myp_sep nan];
    end
end

disp(['thickness=' num2str(nanmean(thick)*1e3) '+/-' num2str(nanstd(thick)*1e3)])
disp(['separation=' num2str(nanmean(myo_myp_sep)*1e3) '+/-' num2str(nanstd(myo_myp_sep)*1e3)])
disp(['bundle fraction=' num2str(nanmean(bndl_act_frac)) '+/-' num2str(nanstd(bndl_act_frac))])
disp(['bundle in myo zone=' num2str(nanmean(bndl_act_myo_frac)) '+/-' num2str(nanstd(bndl_act_myo_frac))])
disp(['bundle in myp zone=' num2str(nanmean(bndl_act_myp_frac)) '+/-' num2str(nanstd(bndl_act_myp_frac))])