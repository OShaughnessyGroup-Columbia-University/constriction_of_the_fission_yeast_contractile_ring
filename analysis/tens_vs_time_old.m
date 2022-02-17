clear
% preallocate
tens_mean_store = [];
tens_sd_store = [];
% go through time
for t = 1:30
    t
    mean_tens_vs_x
    % relative position of 2 micron
    rp2m = 2 / circring;
    % start and end points of averaging
    startpt = round(100 * rp2m);
    endpt = round(100 * (1-rp2m));
    % averaging
    tens_mean_store = [tens_mean_store, mean(mean(tmat(:,[1:startpt,endpt:end])))];
    tens_sd_store = [tens_sd_store, std(mean(tmat(:,[1:startpt,endpt:end]),2))];
end

% plot tension versus time
errorbar(1:length(tens_mean_store)-1, tens_mean_store, tens_sd_store)