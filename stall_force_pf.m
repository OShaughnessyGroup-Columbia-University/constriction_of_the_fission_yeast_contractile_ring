function sfpf = stall_force_pf (bmmat,fone,nsat,nhead,fhead)
% the output sfpf is a full vector.
% will be repeatedly run.

    if size(bmmat,2)~=0  % if there is some myosin
        % calculate stall force per filament of each myosin cluster

        % number of filaments on each myosin cluster
        nfx = sum(bmmat);
        nfx = full(nfx);

        % for Myo2 clusters
        sfpf = (nfx <= nsat) * fone ... if at most 10 filaments on this cluster, stall force is 4pN
            + (nfx > nsat) * nhead * fhead ./ nfx; % if more than 10 filaments on this cluster, a total of 40pN is shared

    %     % for Myp2 clusters
    %     temp = (nfx <= nsat) * fone ... if at most 10 filaments on this cluster, stall force is 4pN
    %         + (nfx > nsat) * nhead * fheadmyp ./ nfx; % if more than 10 filaments on this cluster, a total of 40pN is shared
    %     sfpf = sfpf + (~bancm).*temp;   % Myp2 stall force per filament is fhead.

        % make a column vector
        sfpf = sfpf';
        % replace NaN's by 1e-6
        sfpf(isnan(sfpf)) = 1e-6;
    else % there's no myosin
        sfpf=double.empty();
    end
end