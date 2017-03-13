function [] = manualCorrections(videoPath,refOffset,loessSmooth,minObject,thresh)

    file = './TrackingResuts1A_12-Mar-2017.csv';
    data = csvread(file,1);

    vr = VideoReader(videoPath);
    buffer = 1;
    lastFrame = max(data(:,8));
    %% Options
    bwConnectivity = 8;

    %% Initialize results
    results = [];

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
    startFrame = data(1,8);
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
    imshow(loadBatch(startFrame))
    disp('Select upper left and lower right region where eel is')
    [x,y] = ginput(2); %% User selects tracking box start position
    x = int32(x);
    y = int32(y);
    close(2)
    bbox = [x(1) x(2) y(1) y(2)];
    
    %% Show first frame
    currentFrame = startFrame;
    function showFrame()
        
        idxFrame = find(data(:,8) == currentFrame);
        
        %% Define bbox location
        width = abs(bbox(2) - bbox(1))/2;
        dataTemp = data(idxFrame,9:10);x
        dataTemp = dataTemp(find(~(isnan(dataTemp(:,1)) | isnan(dataTemp(:,2)))),:);
        bbox([1 2]) = [-width width] + mean(dataTemp(:,1));
        height = abs(bbox(3) - bbox(4))/2;
        bbox(3:4) = [-height height] + mean(dataTemp(:,2));
        
        %% Plot results
        % Show detectons (fitted vs pixel centers)
        currentImage = loadBatch(currentFrame);
        figure(1);
        hold off;
        imshow(currentImage(bbox(3):bbox(4),bbox(1):bbox(2))); hold on;
        %% Overlay midline
        scatter(data(idxFrame,9) - double(bbox(1)) + 1,...
            data(idxFrame,10) - double(bbox(3)) + 1,...
            [],[1 0 0],'o'); 
        plot(data(idxFrame,9) - double(bbox(1)) + 1,...
            data(idxFrame,10) - double(bbox(3)) + 1,'color',[1 0 0],'linewidth',2);
        title('Smoothed midline overlayed on source image');
        text(10,10,['Frame: ' num2str(currentFrame)]);
    end

    showFrame();

    %% Set keyboard input
    set(gcf,'KeyPressFcn', @KeyPress);
    function KeyPress(ObjH, EventData)
        Key = get(ObjH, 'CurrentCharacter');
        disp(['Key: '  Key]);
        switch Key
            case 'z'
                frms = unique(data(:,8));
                idx = find(frms  == currentFrame);
                if idx > 1
                    currentFrame = frms(idx - 1);
                    showFrame();
                end
            case 'x'
                frms  = unique(data(:,8));
                idx = find(frms  == currentFrame);
                if idx < length(frms)
                    currentFrame = frms(idx + 1);
                    showFrame();
                end
            case 'a'
                disp('Click on the image twice.')
                disp('First the head, then the tail.')
                disp('Then the midline will be recalculated')
                [x,y] = ginput(2);
                idxFrame = find(data(:,8) == currentFrame);
                data(idxFrame,3:4) = repmat([x(1),y(1)],length(idxFrame),1);
                data(idxFrame,5:6) = repmat([x(2),y(2)],length(idxFrame),1);
                motionDetection;
            case 's'
                disp('Use the mouse to reposition midline points manually.')
                disp('Press enter when finished, then the midlien will be redrawn')
                %% Draw midpoint lines on image
            case 'o'
                disp('Results written to new file!')
        end
    end


    
    %% Loop through all frames
    function motionDetection();
        previousDetections = 0;
        i = currentFrame;
        %% Load current frame
        currentImage = loadBatch(i);
        idxFrame = find(data(:,8) == currentFrame);
        
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
            %% Load detected objects larger than minimum object size
            for j = 1:length(objects.PixelIdxList)
                temp = length(objects.PixelIdxList{j});
                if temp > minObject;
                    bwbb(objects.PixelIdxList{j}) = 1;
                end
            end

            bw(rows,cols) = bwbb;

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
            x = sqrt((transposedPoints(:,1) - detections(:,4)).^2 + (transposedPoints(:,2) - detections(:,5)).^2);
            y = distanceFromMidline;
            % Remove last _ points
            if ignoreLast > 0
                x((end - ignoreLast + 1):end,:) = NaN;
                y((end - ignoreLast + 1):end,:) = NaN;
            end
            curve = [x smooth(x,y,loessSmooth,'rloess')];

            %% Interpolate missing points (cubic spline)
            lastPoint = size(curve,1) - ignoreLast;
            x = ((0:(midlineResolution - ignoreLast))*res);
            temp = interp1(curve(1:lastPoint,1),curve(1:lastPoint,2),...
                x,'spline');
            curve = [x' temp'];

            %% Calc length of smooth midline
            midlineLength = sum(sqrt(sum(diff(curve).^2)));

            %% Transform normalized points back to origional image
            xDisp= curve(:,1) * cos(atan(mMidline));
            yDisp= -curve(:,2) + (curve(:,1) * tan(atan(mMidline)));
            absXSmooth = xDisp + midline(1,1);
            absYSmooth = midline(1,2) + yDisp;

            %% Append NaNs for ignored points
            curve = [curve ;repmat(NaN,ignoreLast,2)];
            absXSmooth = [absXSmooth;repmat(NaN,ignoreLast,1)];
            absYSmooth= [absYSmooth ;repmat(NaN,ignoreLast,1)];

            %% Store results
            points = size(curve,1);
            data(idxFrame,:) = [curve ...  % Relative to swim path
                repmat(midline(1,1),points,1) repmat(midline(1,2),points,1 ),... % Head
                repmat(midline(2,1),points,1) repmat(midline(2,2),points,1 )...  % Tail
                repmat(midlineLength,points,1) ... % Midline Length
                repmat(i,points,1) ... % Frame
                absXSmooth absYSmooth ...  % Absolute midline position
                repmat(midlineResolution,points,1)]; % Midline Points

            %% Convert midline to swimming vector
            % Find first swim position
            idx = find(data(:,8) == startFrame);
            tempX = data(idx,9);
            tempY = data(idx,10);
            idx = find(~(isnan(tempX) | isnan(tempY)));
            start = [mean(tempX(idx)) mean(tempY(idx))];
            % Find last  swim position
            idx = find(data(:,8) == lastFrame);
            tempX = data(idx,9);
            tempY = data(idx,10);
            idx = find(~(isnan(tempX) | isnan(tempY)));
            finish = [mean(tempX(idx)) mean(tempY(idx))];

            % Translocate midline
            dataTemp = data(idxFrame,:);
            y = [];
            for k = 1:(size(idxFrame,1)) 
                y = [y; point2LineDistance(start,finish,dataTemp(k,9:10))];
            end
            % Define swimming path line
            m = (finish(2) - start(2)) / (finish(1) - start(1));
            if m == 0; m = 0.0000001;end;
            if m == -Inf; m = 999999;end;
            b = start(2) - m*start(1);

            [xSnap,ySnap] = snapPointsToLine(dataTemp(:,9),dataTemp(:,10),m,b);

            x = sqrt((xSnap-start(1)).^2 + (ySnap-start(2)).^2);
            % Correct negative values
            idx =  find(start(2) > (m*xSnap + b));
            x(idx) = -x(idx);

            data(idxFrame,9:10) = [x y];
        end
    end
 
    
    %% Save results
    function saveResults
        fw = fopen([resultsFile '_corrected.csv'],'w');
        fprintf(fw,'%s \n','CurveX,CurveY,HeadX,HeadY,TailX,TailY,MidlineLength,Frame,AbsX,AbsY,MidlineSlope');
        fprintf(fw,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n',data');
        fclose(fw);
    end
    
end

