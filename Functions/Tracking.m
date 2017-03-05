function Tracking(videoPath,firstDetection,lastFrame,refOffset,...
    midlineResolution,loessSmooth,buffer,frameInterval,...
    thresh,Individual, Treatment,video)
    vr = VideoReader(videoPath);
    addpath('Functions')
    
    %% Options
    bwConnectivity = 8;

    %% Initialize results
    results = [];

    if video == 1;
        vw = VideoWriter(['TrackingVideo' Individual Treatment '_' date '.avi']);
        vw.FrameRate = 5;
        open(vw);
    end
    %% Function for loading video frames (current frame)
    % This function facilitates loeading multiple frames at once.
    % For computers with more memory, loading frames in batches will speed
    % up processing time on video formats other than MJPEG
    function [frame] = loadBatch(idx)
      
        %% Check that idx is within video boundaries
        if idx > lastFrame || idx < 1
            % Return a white frame (full detection) if frame number is out
            % of bounds
           frame = ones(vr.Height,vr.Width);
           frame(:,:) = 255;
           frame = uint8(frame);
        else
            %% Check that new frame is within buffer
            % if not, refresh buffer
            if idx > endFrame || idx < startFrame
                startFrame = max(idx,1);
                endFrame = min([(startFrame + buffer - 1) vr.NumberOfFrames]);
                frames = read(vr,[startFrame,endFrame]);
            end
            % Load frames and convert to graycale
            frame = rgb2gray(frames(:,:,:,idx - startFrame + 1));
        end
    end
    
    %% Load reference frame
    startFrame = 1;
    endFrame = startFrame + buffer - 1;
    frames = read(vr,[startFrame,endFrame]);
    % The reference frame is redefined on each frame of video.  It is the
    % region around the eel n frames before or after the current frame,
    % where n = refOffset.  This is useful in case the camera moves
    % slightly during the setup, then only a portion of the video has
    % tracking errors.  This can be improved so that both the reference
    % before and after the current frame is used, and the one with the
    % least motion detection is used.

    %% First frame where eel is fully visible
    figure(2);
    imshow(loadBatch(firstDetection))
    disp('Select upper left and lower right region where eel is')
    [x,y] = ginput(2); %% User selects tracking box start position
    x = int32(x);
    y = int32(y);
    close(2)
    bbox = [x(1) x(2) y(1) y(2)];
    
    %% Loop through all frames
    previousDetections = 0;
    framesToAnalyse = firstDetection:frameInterval:min(((vr.NumberOfFrames) -1),lastFrame)
    for i= framesToAnalyse
        %% Load current frame
        currentImage = loadBatch(i);
       
        %% Load reference image (both before or after current frame)
        refImageBefore = loadBatch(i - refOffset);
        refImageAfter = loadBatch(i + refOffset);
        
        %% Motion detection
        [bwRaw,greyScale] = DetectFish(refImageBefore,refImageAfter,...
            currentImage,bbox,thresh);
        
        %% Image Cleaning
        % remove long thin detected areas
        bw = bwmorph(bwRaw, 'majority');
        
        %% Grab largest continuous object from motion
        rows = bbox(3):bbox(4);
        cols = bbox(1):bbox(2);
        objects = bwconncomp(bw(rows,cols),bwConnectivity);
        objects.sizes = cellfun('size', objects.PixelIdxList, 1);
        [~,I] = max(objects.sizes);
        bwbb = bw(rows,cols);
        bwbb(:,:) = 0;
        
        % If no motion detections, skip this frame
        if ~isempty(objects.PixelIdxList);
            % Load largest detected object
             bwbb(objects.PixelIdxList{I}) = 1;
            %% Load detected objects which overlap with previous objects
            % Detected objects which overlap with objects detected in the
            % previous frames will be included here
            for j = 1:length(objects.PixelIdxList)
                temp = intersect(previousDetections,objects.PixelIdxList{j});
                if ~isempty(temp)
                    bwbb(objects.PixelIdxList{j}) = 1;
                end
            end
            bw(rows,cols) = bwbb;
            
            %% Estimate head location
            % Average vector per frame (forward movement)
            
            % 
            
            previousDetections = find(bwbb);
            
            %% Draw Convex hull and calculate midline
            cv = cHull(bw);
            [midline,D] = longestLine(cv);
            if(midline(1,1) > midline(2,1))
                %% Make sure first point is on left side
                midline = midline([2 1],:);
            end
           
            %% Estimate midline points
            [detections,transposedPoints,distanceFromMidline,...
                mMidline,bMidline] = loopMidlineSegments(...
                midline(1,:),midline(2,:),midlineResolution,D,bw);
            
            res = (D/midlineResolution);
            
            %% Add fitted spline (smooth using loess)            
            curve = smooth((detections(:,3)-1)*res,...
                distanceFromMidline,loessSmooth,'rloess');
            curve = [(detections(:,3)-1)*res curve];

            %% Interpolate missing points (cubic spline)
            temp = interp1(curve(:,1),curve(:,2),...
                (0:midlineResolution)*res,'spline');                    
            curve = [((0:midlineResolution)*res)' temp'];
            %% Calc length of smooth midline
            midlineLength = sum(sqrt(sum(diff(curve).^2)));

            %% Transform normalized points back to origional image
            xDisp= curve(:,1) * cos(atan(mMidline));
            yDisp= -curve(:,2) + (curve(:,1) * tan(atan(mMidline)));
            absXSmooth = xDisp + midline(1,1);
            absYSmooth = midline(1,2) + yDisp;
            
            %% Store results
            points = size(curve,1);
            data = [curve ...
                repmat(midline(1,1),points,1) repmat(midline(1,2),points,1 ),...
                repmat(midline(2,1),points,1) repmat(midline(2,2),points,1 )...
                repmat(midlineLength,points,1) repmat(i,points,1) absXSmooth absYSmooth...
                repmat(mMidline,points,1)];
            results = [results; data];
      
            %% Redefine rectangle area (recenter)
            % This updates the motion detection area for the following frame
            c = [((midline(1,1) + midline(2,1))/2) ((midline(1,2) + midline(2,2))/2)];
            bbox = [(c(1) - ((bbox(2) - bbox(1))/2)) ...
                (c(1) + ((bbox(2) - bbox(1))/2)) ...
                c(2) - ((bbox(4) - bbox(3))/2) ...
                c(2) + ((bbox(4) - bbox(3))/2)];

            %% Plot results
            % Show detectons (fitted vs pixel centers)
            subplot(4,1,1);hold off;
            imshow(bw(bbox(3):bbox(4),bbox(1):bbox(2))*0.5 + ...
                bwRaw(bbox(3):bbox(4),bbox(1):bbox(2))*0.5); hold on;
            scatter(detections(:,1) - double(bbox(1)) + 1,...
                detections(:,2) - double(bbox(3)) + 1 ...
                ,[],[0.5 0 0],'o');
            scatter(transposedPoints(:,1) - double(bbox(1)) + 1 ...
                ,transposedPoints(:,2) - double(bbox(3)) + 1 ...
                ,[],[1 0 0],'x');
            title('Midline Pixel Detection');
            % Show normalized and fitted curve
            subplot(4,1,2); hold off;
            scatter((detections(:,3)-1)*res,distanceFromMidline,...
                [],[1 0 0],'filled','o'); 
            xlim([min(data(:,1)) max(data(:,1))]); axis equal;
            hold on; plot(data(:,1),data(:,2)); grid on;
            legend('Average Midline','Smoothed Midline');
            title('Normalized & Smoothed midline')
            % Show image subtraction
            figure(1);subplot(4,1,3);hold off;
            imshow(greyScale(bbox(3):bbox(4),bbox(1):bbox(2)) - currentImage(bbox(3):bbox(4),bbox(1):bbox(2)))            
            title('Image subtraction (Ref image before current frame)'); hold on;
            text(10,10,...
                ['Pixels Detected: ' num2str(length(find(bw)))...
                ' Frame: ' num2str(i)...
                ' Straight Midline Length: ' num2str(D)],...
                'color',[0.8 0 0]);
            % Show origional image overlayed with smoothing spline
            subplot(4,1,4); hold off;
            imshow(currentImage(bbox(3):bbox(4),bbox(1):bbox(2))); hold on;
            scatter(transposedPoints(:,1) - double(bbox(1)) + 1,...
                transposedPoints(:,2) - double(bbox(3)) + 1,...
                [],[1 0 0],'o'); 
            plot(data(:,9) - double(bbox(1)) + 1,...
                data(:,10) - double(bbox(3)) + 1,'color',[1 0 0],'linewidth',2);
            title('Smoothed midline overlayed on source image');
            
            if video == 1;
                 writeVideo(vw,getframe(1));
            end
        end
    end
    
    fw = fopen(['TrackingResuts' Individual Treatment '_' date '.csv'],'w');
    fprintf(fw,'%s \n','CurveX,CurveY,Midline1X,Midline1Y,Midline2X,Midline2Y,SmoothedMidlineLength,Frame,AbsX,AbsY,MidlineSlope');
    fprintf(fw,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n',results');
    fclose(fw);
    if video == 1;
        close(vw);
    end
    disp('Analysis complete!')
 
end