%% decide what files are available
clear
a = dir;
matname = 'ds_4_';
writeObj = VideoWriter('intf_ring4_simb','MPEG-4');
writeObj.Quality = 100;
writeObj.FrameRate = 1;
open(writeObj);
% load file
figure
for iiiload = 0:1e5
    file_name = strcat(matname,num2str(iiiload),'s.mat')
    if exist(file_name, 'file')
        load(strcat(matname,num2str(iiiload),'s.mat'), 'rmyo', 'rbead', 'ifor', 'ipt', 'bancf', 'bancm','nbead','xmat')
        ringplot_det_ring
        axis([-.3 .8 1.8 2.2]);
        %axis([-2.5, 2.5, -2.5, 2.5]);
        title(strcat('t = ',num2str(iiiload),' s'))
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        hold on
        %viscircles([0,0], 2.13,'EdgeColor','c','LineWidth',6);
        F(iiiload+1) = getframe(gcf);
        writeVideo(writeObj, F(iiiload+1));
        close
    else 
        break
    end
end
close(writeObj);