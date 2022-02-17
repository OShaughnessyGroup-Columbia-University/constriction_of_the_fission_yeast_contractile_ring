function [bndl_mean, bndl_std, thick_mean, thick_std, width_mean,...
    width_std, sep_mean, sep_std] = rHistogramTime(tag, tlist, unit)

bndl_act_frac = [];
bndl_act_myo_frac = [];
bndl_act_myp_frac = [];
bndl_myp_frac = [];
thick = [];
width = [];
rbndl = [];
xsect_fil = [];
myo_myp_sep = [];

% tlist = 0:10:tmax;
% unit = 'sec';

% tlist = 0:tmax;
% unit = 'min';
for it = 1:length(tlist)
    snap_time = tlist(it)
    for k=1:64
        if exist(sprintf(['%s/%s_%d_' num2str(snap_time) unit '.mat'],tag,tag,k)) == 2
            load(sprintf(['%s/%s_%d_' num2str(snap_time) unit '.mat'],tag,tag,k))
            rHist_noplot
            %break;
            bndl_act_frac = [bndl_act_frac cnt/size(rbead,2)*100];
            bndl_act_myo_frac = [bndl_act_myo_frac cntMyo/cnt*100];
            bndl_act_myp_frac = [bndl_act_myp_frac cntMyp/cnt*100];
            bndl_myp_frac = [bndl_myp_frac cntMypInBundle/size(rmyp,2)*100];

            % 95 thickness
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
           disp(['No file i=' num2str(k) ' t=' num2str(tlist(it))])
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
    width_mean(it) = nanmean(width);
    width_std(it) = nanstd(width);
    thick_mean(it) = nanmean(thick);
    thick_std(it) = nanstd(thick);
    sep_mean(it) = nanmean(myo_myp_sep);
    sep_std(it) = nanstd(myo_myp_sep);
    bndl_mean(it) = nanmean(bndl_act_frac);
    bndl_std(it) = nanstd(bndl_act_frac);
end
end

% disp('********************************************************************************')
% disp(['thickness=' num2str(nanmean(thick)) '+/-' num2str(nanstd(thick))])
% disp(['separation=' num2str(nanmean(myo_myp_sep)) '+/-' num2str(nanstd(myo_myp_sep))])
% disp(['bundle fraction=' num2str(nanmean(bndl_act_frac)) '+/-' num2str(nanstd(bndl_act_frac))])
% disp(['bundle in myo zone=' num2str(nanmean(bndl_act_myo_frac)) '+/-' num2str(nanstd(bndl_act_myo_frac))])
% disp(['bundle in myp zone=' num2str(nanmean(bndl_act_myp_frac)) '+/-' num2str(nanstd(bndl_act_myp_frac))])