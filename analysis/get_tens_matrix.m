tens_cell = {};
for iload = 1:64
    tens_vec = [];
    for jload = 1:30
        filename = ['sizep5_', num2str(iload), '_', num2str(jload),'min.mat'];
        if exist(filename,'file')
            load(filename)
%             tension_laplace_v1
            tension_circ
            tens_vec = [tens_vec, mean(tens)];
        else
            continue
        end
    end
    tens_cell{iload} = tens_vec;
end
save('temp.mat','tens_cell')