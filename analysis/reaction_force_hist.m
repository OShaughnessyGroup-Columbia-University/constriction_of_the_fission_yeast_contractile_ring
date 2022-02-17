figure
hold on
for irun = 18:19
    fp_phi = [];
    try
    load(['noagm_' num2str(irun) '_00sec.mat'])
    for it = 30:30:600
        try
        load(['noagm_' num2str(irun) '_' num2str(it) 'sec.mat'])
        force_enum
        close
        close
        fp_phi = [fp_phi; fpull(:, 2)];
        catch
            continue
        end
    end
    histogram(fp_phi, -30:30, 'normalization', 'probability', 'displayname',...
        ['irun = ' num2str(irun)])
    catch
        continue
    end
end