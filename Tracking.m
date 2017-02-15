function Tracking(videoPath,firstDetection)
%   startFrame = 150
%   video = 'ExampleVideo/sprint_20160722_024359_PM_20160722_025355_PM.avi';
    vr = VideoReader(videoPath);
    addpath('Functions')
    
    %% Options
    bwConnectivity = 8;
    buffer = 600;  % load multiple video frames at a time (increases processing speed)
    % How many frames ahead (or behind) the reference image will be grabbed from
    refOffset = 150;
    minSize = 50; % Minimum size for a detected object
    midlineResolution = 10;
    curvePoints = 10;
    %% Initialize results
    results = [];

    %% Function for loading video frames (current frame)
    % This function facilitates loeading multiple frames at once.
    % For computers with more memory, loading frames in batches will speed
    % up processing time on video formats other than MJPEG
    function [frame] = loadBatch(idx)
        %Returns an int8 matrix of greyscale values
        if idx > vr.NumberOfFrames || idx < 1
            % Return a white frame (full detection) if frame number is out
            % of bounds
           frame = ones(vr.Height,vr.Width);
           frame(:,:) = 255;
           frame = uint8(frame);
        else
            if idx > endFrame || idx < startFrame
                startFrame = idx;
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
    

    
    %% loop through all frames
    for i=firstDetection:((vr.NumberOfFrames) -1)
        %% Load current frame
        currentImage = loadBatch(i);
       
        %% Load reference image (both before or after current frame)
        refImageBefore = loadBatch(i - refOffset);
        refImageAfter = loadBatch(i + refOffset);
        
        %% Create binary image
        % Do background subtraction from both before and after current
        % frame.  Pick the resulting image with the least detected motion.
        % This should avoid detection issues relatied to camera movement
        % and other moving objects eerarlier and later in the video.
        bwBefore = backgroundSubtraction(refImageBefore,currentImage, bbox);
        bwAfter = backgroundSubtraction(refImageAfter,currentImage, bbox);
        if length(find(bwBefore == 1)) <= length(find(bwAfter == 1))
            bw = bwBefore;
        else
            bw = bwAfter;
        end
        
        %% Image Cleaning
        % remove long thin detected areas
        bw = bwmorph(bw, 'majority');
        
        %% Grab largest continous object from motion
        % If largest object is smaller than minSize, use the two largest
        % objects.  This should solve issues related to the eel swimming
        % past a break in the mirrors or tank
        objects = bwconncomp(bw,bwConnectivity);
        objects.sizes = cellfun('size', objects.PixelIdxList, 1);
        [~,I] = sort(objects.sizes,2,'descend');
        bw(:,:) = 0; % Clear image
        % redraw image with only largest object
        bw(objects.PixelIdxList{I(1)}) = 1;
        % Adds second largest object if detection is too small
        if (objects.sizes(I(1,1)) < minSize) && (length(I) > 1)
            bw(objects.PixelIdxList{I(2)}) = 1;
        end
        
        %% Draw Convex hull
        cv = bwboundaries(bw);
        CH = convhull(cv{1}(:,1),cv{1}(:,2));
        cv = cv{1}(CH,[2 1]);
        %% Get midline
        [midline,D] = longestLine(cv);
        
        %% Get detected pixels perpendicular to midline
        disp('Midline detections')
        n = floor(D/midlineResolution); %Get numper of points to plot
        res = D/n; % interval between perpendicular lines
        % Get slope ( y1 - y2 / x1 - x2 )
        mMidline = (midline(2,2) - midline(1,2)) / (midline(2,1) - midline(1,1));
        m = -1/mMidline; % Perpendicular slope
        bStart = midline(1,2) - m*midline(1,1); %intercept of first perpendicular line
        dx = res/(1+(mMidline^2));
        dy = m*dx;
        
        % get coords of detected pixels
        [row,col,~] = find(bw);
        % Plot detected areas 
        figure(1);hold off;imshow(bw);axis equal;hold on;
        plot(midline(:,1),midline(:,2),'color',[0 1 0]);
        xlim(bbox(1:2)); ylim(bbox(3:4));
        % Loop through all possible line segments
        detections = [];
        figure(2)

        for j=0:n
            % plot line
            b = bStart - (j*dy) - (1/dx);
            plot(bbox(1:2),m*bbox(1:2) + b,'color',[0 1 0])
            % Find intersecting pixels
            x = [(col - 0.5) (col + 0.5)];
            y = [(row - 0.5) (row + 0.5)];
            % predict y values for intersecting points
            yPred = m*x + b;
            % Save intersected pixels
            [idx,~,~] = find((y(:,1) <= max(yPred,[],2)) & (y(:,2) >= min(yPred,[],2)));
            if size(idx,1) > 0
                detections = [detections; col(idx) row(idx) repmat(j+1,length(idx),1)];
                % move detections onto line
            end
        end
        
        % Normalize image (center around 0,0)
        detectionsRange = [min(detections(:,1)) max(detections(:,1)) min(detections(:,2)) max(detections(:,2))];
        detectionsNorm = [(detections(:,1) - mean(midline(:,1))) (detections(:,2) - mean(midline(:,2)))];
        % Rotate data to be horizontal
        
        %% Show detectons
        scatter(detections(:,1),detections(:,2),[],[1 0 0],'filled');
        
        
        %% Redefine rectangle area (recenter)
        %% This updates the motion detection area for the following frame
        c = [((midline(1,1) + midline(2,1))/2) ((midline(1,2) + midline(2,2))/2)];
        bbox = [(c(1) - ((bbox(2) - bbox(1))/2)) ...
            (c(1) + ((bbox(2) - bbox(1))/2)) ...
            c(2) - ((bbox(4) - bbox(3))/2) ...
            c(2) + ((bbox(4) - bbox(3))/2)];
        
        %% Fit smooth line
        % locations of pixels
        try
            temp = bw2(bbox(3):bbox(4),bbox(1):bbox(2));
        catch
            disp('Cannot load bounding box, eel is likely on the edge of the screen.  Ending analysis')
            break
        end
        idx = find(temp);
        [j,k] = ind2sub(size(temp),idx);
        % Fit line
        [x,y] = fitspline(j,k);
        
        results = [results; [ones(size(x,2),1)*i (int32(x) + bbox(1))' (int32(y) + bbox(3))']];
        
        figure(4)
        subplot(5,1,1);
            imshow(abs(refImageBefore(bbox(3):bbox(4),bbox(1):bbox(2)) -...
                currentImage(bbox(3):bbox(4),bbox(1):bbox(2))))            
            title('Image subtraction (Ref image before current frame)')
        subplot(5,1,2);
            imshow(abs(refImageBefore(bbox(3):bbox(4),bbox(1):bbox(2)) -...
                currentImage(bbox(3):bbox(4),bbox(1):bbox(2))))            
            title('Image subtraction (Ref image after current frame)')
        subplot(5,1,3);
            imshow(bw(bbox(3):bbox(4),bbox(1):bbox(2)));
            title('Small objects removed')
        subplot(5,1,4);
            imshow(currentImage(bbox(3):bbox(4),bbox(1):bbox(2))); hold on
            % transparent color mask
            mask = ones(size(temp));
            mask(:,:,2) = ones(size(temp))*0.5;
            mask(:,:,3) = ones(size(temp))*0.5;
            hmask = imshow(mask);
            set(hmask,'AlphaData',temp*0.5)
            hold off;
            title('Cleaned object overlayed on video')
        subplot(5,1,5);
            imshow(currentImage(bbox(3):bbox(4),bbox(1):bbox(2))); hold on;
            scatter(x,y,'.r','lineWidth',2); 
            text(1,1,['Frame: ' num2str(i) '  Seconds: ' num2str(i/vr.FrameRate)],'color','red');hold off
            title('Averaged y-values (estimate of center line)')
    end
    
    csvwrite(['TrackingResuts' date '.csv'], results)
    disp('Analysis complete!')
end