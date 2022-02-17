function rdiff = get_rdiff(rbead,~,ipt)
% will be run repeatedly

    % vector from this bead to the next bead
    rdiff = circshift(rbead,[0,-1]) - rbead;
    % normalize 
    for i = 1:size(rdiff,2)
        rdiff(:,i) = rdiff(:,i) / norm(rdiff(:,i)); 
    end
    % do not count rdiff on ipt
    rdiff(:,ipt) = zeros(3,numel(ipt));
    rdiff = -rdiff; % for pted end anchoring
end
