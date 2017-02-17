function writeEelPreview(file, fps)
    vw = VideoWriter('EelPreview.avi');
    vw.FrameRate = fps;
    data = csvread(file,1,0);
    
    xmin = min(data(:,9));
    xmax = max(data(:,9));
	ymin = min(data(:,10));
    ymax = max(data(:,10));
    
    frames = unique(data(:,8));
    
    %% Intialize plot
    hf = figure(1);
    hscatter = scatter(0,0,'r.');
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    axis equal
    
    open(vw)
    %% Plot poinits
    for i = frames' 
        idx = find(data(:,8) == i);
        set(hscatter,'XData',data(idx,9))
        set(hscatter,'YData',data(idx,10))
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        
        writeVideo(vw,getframe(1))
    end
    close(vw)
    beep;
    disp('Video Written!  :P');
end