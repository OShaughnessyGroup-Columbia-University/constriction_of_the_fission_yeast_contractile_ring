clear
% preallocate
ilist = 1:10;
tlist = 1:35;
tens_store = nan(numel(ilist), numel(tlist));
for itemp = 1:numel(ilist)
    for t = tlist
        filename = strcat('fcnew_', num2str(ilist(itemp)), '_', num2str(t), 'min.mat');
        if exist(filename,'file')
                load(filename)
            tension_circ
            % averaging
            tens_store(itemp,t) = mean(tens);
        else
            break
        end
    end
end
bad_cols = any(isnan(tens_store));
tens_store(:,bad_cols) = [];
tlist(bad_cols) = [];

tens_store = tens_store * 1.5;

tens_mean = mean(tens_store);
tens_sd = std(tens_store);
%% plot tension versus phi
r_ring_t = 1.85 - 0.07*(tlist-10);
phi = acosd(r_ring_t / 1.85);
phi(real(phi)==0)=nan;
errorbar(real(phi),tens_mean,tens_sd,'-ok')
xlabel('\phi (deg)','FontSize',18)
ylabel('Tension (pN)','FontSize',18)
% axis([0,26,0,1600])