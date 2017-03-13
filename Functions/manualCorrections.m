function [] = manualCorrections(videoPath,refOffset,loessSmooth,minObject,thresh)

    [file,path] = uigetfile({'*.*','All Files'},'Select Results File');
    resultsFile = [path,file];
    data = csvread(resultsFile,1);

    vr = VideoReader(videoPath);
    buffer = 1;
    lastFrame = max(data(:,8));
    %% Options
    bwConnectivity = 8;

    %% Initialize results
    midlineResolution = data(1,11);
    ignoreLast = data(1,12);

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
    firstDetection = startFrame;
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
        dataTemp = data(idxFrame,9:10);
        dataTemp = dataTemp(find(~(isnan(dataTemp(:,1)) | isnan(dataTemp(:,2)))),:);
        bbox([1 2]) = [-width width] + mean(dataTemp(:,1));
        height = abs(bbox(3) - bbox(4))/2;
        bbox(3:4) = [-height height] + mean(dataTemp(:,2));
        
        %% Plot results
        % Show detectons (fitted vs pixel centers)
        currentImage = loadBatch(currentFrame);
        figure(3);
        hold off;
        imshow(currentImage(bbox(3):bbox(4),bbox(1):bbox(2))); hold on;
        %% Overlay midline
        scatter(data(idxFrame,9) - double(bbox(1)) + 1,...
            data(idxFrame,10) - double(bbox(3)) + 1,...
            [],[1 0 0],'o'); 
        plot(data(idxFrame,9) - double(bbox(1)) + 1,...
            data(idxFrame,10) - double(bbox(3)) + 1,'color',[1 0 0],'linewidth',2);
        scatter(data(idxFrame(1),[3 5]) - double(bbox(1)), data(idxFrame(1),[4 6]) - double(bbox(3)),...
            [],[1 0 0],'og'); hold off;
        
        title('Smoothed midline overlayed on source image');
        text(10,10,['Frame: ' num2str(currentFrame)]);
        set(gcf,'KeyPressFcn', @KeyPress);
    end

    showFrame();

    %% Set keyboard input
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
                data(idxFrame,3:4) = repmat([x(1) + double(bbox(1)),y(1) + double(bbox(3))],length(idxFrame),1);
                data(idxFrame,5:6) = repmat([x(2) + double(bbox(1)),y(2) + double(bbox(3))],length(idxFrame),1);
                motionDetection;
            case 's'
                disp('Use the mouse to reposition midline points manually.')
                disp('Press enter when finished, then the midlien will be redrawn')
                %% Draw midpoint lines on image
                manualDetection;
            case 'o'
                saveResults;
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
            %% Load detected objects larger than minimum object size or overlap head or tail point
            for j = 1:length(objects.PixelIdxList)
                temp = length(objects.PixelIdxList{j});
                if temp > minObject;
                    bwbb(objects.PixelIdxList{j}) = 1;
                end
            end

            bw(rows,cols) = bwbb;

            %% Draw Convex hull and calculate midline
            midline = [data(idxFrame(1),[3 4]);data(idxFrame(1),[5 6])];
            D = sqrt((midline(1,1) - midline(2,1)).^2 + (midline(1,2) - midline(2,2)).^2);
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
                repmat(midlineResolution,points,1) ... % midlien resolution
                repmat(ignoreLast,points,1)]; % Ignore last n points

            %% Convert midline to swimming vector
            % Find first swim position
            idx = find(data(:,8) == firstDetection);
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

            data(idxFrame,1:2) = [x y];
            showFrame() %% Refresh image
        end
    end

%% Loop through all frames
    function manualDetection();
        
        %% Show green lines for snapping on image
        % Get head and tail positions
        idxFrame = find(data(:,8) == currentFrame);
        midline(1,:) = data(idxFrame(1),3:4);
        midline(2,:) = data(idxFrame(1),5:6);
        D = sqrt((midline(1,1) - midline(2,1)).^2 + (midline(1,2) - midline(2,2)).^2);
        
        % Define lines
        [mOut,bOut] = DefineMidlineSegments(midline(1,:) - double([bbox(1) bbox(3)]) + 1,midline(2,:) - double([bbox(1) bbox(3)]) + 1,midlineResolution,D);
        figure(3)
        Q1 = [];
        Q2 = [];
        for i = 1:length(mOut)
            hold on;
            xLine = [0;bbox(2) - bbox(1)];
            yLine = [mOut(i)*(xLine(1)) + bOut(i);...
                mOut(i)*(xLine(2)) + bOut(i)];
            plot(xLine,yLine,'g');hold off;
            Q1(i,:) = [xLine(1) yLine(1)];
            Q2(i,:) = [xLine(2) yLine(2)];
        end

        
        %% Replace points on closest line (get user input)
        disp('Click on the green lines where youd like to correct the points');
        disp('When finished, press enter')
        disp('After, your clicked lcoations will snap to the lines and become the new locations')
        [xIn,yIn] = ginput();
        
        % Curve is absolute position
        curve = data(idxFrame,9:10);
        % convert to bounding box location
        curve = [(curve(:,1) - double(bbox(1)) + 1 ) (curve(:,2) + 1 - double(bbox(3)))]; 
        % scatter(curve(:,1),curve(:,2))
        
        % Find closest line to each point
        manualPoints = [];
        for i = 1:length(xIn)
            distances = [];
            for j = 1:(length(mOut) - ignoreLast)
                distances(j) = point2LineDistance(Q1(j,:),Q2(j,:),[xIn(i) yIn(i)]);
            end
            % Get closest line
            [~,I] = min(abs(distances));
            manualPoints = [manualPoints I];
            % Snap to closest line
            [xSnap,ySnap] = snapPointsToLine(xIn(i),yIn(i),mOut(I),bOut(I));
            curve(I,:) = [xSnap,ySnap];
        end
        manualPoints = unique(manualPoints);
        % scatter(curve(:,1),curve(:,2))
        
        % Convert back to absolute position
        curve = [(curve(:,1) - 1 + double(bbox(1))) (curve(:,2) - 1 + double(bbox(3)))]; 
        
        %% Translocate to major axis reference
        y = [];
        for k = 1:(size(idxFrame,1)) 
            y = [y; point2LineDistance(midline(1,:),midline(2,:),curve(k,:))];
        end
        % Define swimming path line
        m = (midline(2,2) - midline(1,2)) / (midline(2,1) - midline(1,1));
        if m == 0; m = 0.0000001;end;
        if m == -Inf; m = 999999;end;
        b = midline(1,2) - m*midline(1,1);

        [xSnap,ySnap] = snapPointsToLine(curve(:,1),curve(:,2),m,b);
        %scatter(xSnap,ySnap)

        x = sqrt((xSnap-midline(1,1)).^2 + (ySnap-midline(1,2)).^2);
        % Correct negative values
%        idx =  find(midline(1,2) > (m*xSnap + b));
        % x(idx) = -x(idx);
        % scatter(x,y)
        curve = [x y];
        
        %% Interpolate range of manually specified points
        queryPoints = min(manualPoints):max(manualPoints);
        curve(queryPoints,2) = interp1(curve(manualPoints,1),curve(manualPoints,2),curve(queryPoints,1),'spline');

        %% Calc length of smooth midline
        midlineLength = sum(sqrt(sum(diff(curve).^2)));

        %% Transform smoothed and interpolated points back to origional image
        mMidline = (midline(2,2) - midline(1,2)) / (midline(2,1) - midline(1,1));
        xDisp= curve(:,1) * cos(atan(mMidline));
        yDisp= -curve(:,2) + (curve(:,1) * tan(atan(mMidline)));
        absXSmooth = xDisp + midline(1,1);
        absYSmooth = midline(1,2) + yDisp;
        % scatter(absXSmooth,absYSmooth)
        
%        curve = [curve; repmat(nan,ignoreLast,2)];
%        absXSmooth = [absXSmooth; repmat(nan,ignoreLast,1)];
%        absYSmooth= [absYSmooth; repmat(nan,ignoreLast,1)];
        %% Store results
        points = size(curve,1);
        data(idxFrame,:) = [curve...  % Relative to swim path
            repmat(midline(1,1),points,1) repmat(midline(1,2),points,1 ),... % Head
            repmat(midline(2,1),points,1) repmat(midline(2,2),points,1 )...  % Tail
            repmat(midlineLength,points,1) ... % Midline Length
            repmat(currentFrame,points,1) ... % Frame
            absXSmooth  absYSmooth  ...  % Absolute midline position
            repmat(midlineResolution,points,1) ... % midlien resolution
            repmat(ignoreLast,points,1)]; % Ignore last n points

        %% Convert midline to swimming vector
        % Find first swim position
        idx = find(data(:,8) == firstDetection);
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

        data(idxFrame,1:2) = [x y];
        showFrame() %% Refresh image
    end
 
    
    %% Save results
    function saveResults
        fw = fopen([resultsFile '_corrected.csv'],'w');
        fprintf(fw,'%s \n','CurveX,CurveY,HeadX,HeadY,TailX,TailY,MidlineLength,Frame,AbsX,AbsY,MidlinePoints,IgnoreLastNPoints');
        fprintf(fw,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n',data');
        fclose(fw);
    end
    
end

