function writeEelPreview(file,width, height, fps)
    vw = VideoWriter('EelPreview.avi');
    vw.FrameRate = fps;
    data = csvread(file);
    
    xmin = min(data(:,2));
    xmax = max(data(:,2));
	ymin = min(data(:,3));
    ymax = max(data(:,3));
    
    frameStart = min(data(:,1));
    frameEnd = max(data(:,1));
    
    %% Intialize plot
    hf = figure(1);
    hscatter = scatter(0,0,'r.');
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    get(hf)
    set(hf,'Position',[0,0,width,height])
    axis equal
    
    open(vw)
    %% Plot poinits
    for i = frameStart:frameEnd
        idx = find(data(:,1) == i);
        set(hscatter,'XData',data(idx,2))
        set(hscatter,'YData',data(idx,3))
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        
        writeVideo(vw,getframe(1))
    end
    close(vw)
    beep;
    disp('Video Written!  :P');
end