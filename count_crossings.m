function ncross_seg = count_crossings(pfx, irun, it)

ncross_seg = crossDetect(pfx, irun, it);
save([ pfx '_ncross_' num2str(irun) '.mat'])
disp(' * * * * * * * * * * Done * * * * * * * * * * ')

end
