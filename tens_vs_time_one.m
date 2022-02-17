clear
% preallocate
tens_store = [];
nbead_store = [];
% go through time
for t = 1:1e5
    filename = strcat('poob_38_', num2str(t), 'min.mat');
    if exist(filename,'file')
            load(filename)
%         tension_circ
        tension_laplace_v2
        % averaging
        tens_store = [tens_store, mean(tens)];
        nbead_store = [nbead_store, size(rbead,2)];
    else
        break
    end
end
t = 1:t-1;
% plot tension versus time
plot(t,tens_store)
hold on
%plot(t,nbead_store-mean(nbead_store) + mean(tens_store))
hold off