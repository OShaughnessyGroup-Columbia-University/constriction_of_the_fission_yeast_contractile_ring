% calculate the force exerted on the membrane by anchors

% change fc into a 3-by-#seg matrix
fcmat = vec2mat(fc,3)';

% preallocate
segtens = zeros(3,nbead); 

% go thru every bead on the filament
for j = nbead:-1:1
    % if at a pointed end, use fc directly
    if ismember(j,ipt)
        segtens(:,j) = fcmat(:,j);
    else
        segtens(:,j) = fcmat(:,j) + segtens(:,j+1); % segtens on formins are membrane forces
    end
end

% concatenate formins and myosin
membforce = [segtens(:,ifor), fcmat(:,nbead+1:end)];

sum(sqrt(sum(membforce .* membforce))) / 2 / pi