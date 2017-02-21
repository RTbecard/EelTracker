function Tracking(videoPath,firstDetection,lastFrame,refOffset,...
    midlineResolution,loessSmooth,maxMidlineLength,buffer,frameInterval,...
    thresh,Individual, Treatment,video,minChunkSize)
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
        %Returns an int8 matrix of greyscale values
        if idx > endFrame || idx < startFrame
            % Return a white frame (full detection) if frame number is out
            % of bounds
           frame = ones(vr.Height,vr.Width);
           frame(:,:) = 255;
           frame = uint8(frame);
        else
            if idx > lastFrame || idx < startFrame
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
    
    %% loop through all frames
    for i=firstDetection:frameInterval:min(((vr.NumberOfFrames) -1),lastFrame )
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
        bwBefore = backgroundSubtraction(refImageBefore,currentImage, bbox,thresh);
        bwAfter = backgroundSubtraction(refImageAfter,currentImage, bbox,thresh);
        if length(find(bwBefore == 1)) <= length(find(bwAfter == 1))
            bw = bwBefore;
            greyScale = refImageBefore; 
        else
            bw = bwAfter;
            greyScale = refImageAfter; 
        end
        
        bwSubtract = bw;
        
        %% Image Cleaning
        % remove long thin detected areas
        bw = bwmorph(bw, 'majority');
        
        %% Grab largest continous object from motion
        % If largest object is smaller than minSize, use the two largest
        % objects.  This should solve issues related to the eel swimming
        % past a break in the mirrors or tank
        objects = bwconncomp(bw,bwConnectivity);
        objects.sizes = cellfun('size', objects.PixelIdxList, 1);
        
        if ~isempty(objects.PixelIdxList);
            % If no motion detections, skip
            
            [~,I] = sort(objects.sizes,2,'descend');
            bw(:,:) = 0; % Clear image
            % redraw image with only largest object
            bw(objects.PixelIdxList{I(1)}) = 1;
            % Calc midline length
            cv = cHull(bw);
            [midline,D] = longestLine(cv);
            % Adds second largest object if detection is too small
            if ((length(I) > 1) && (D < maxMidlineLength))
                temp = bw;
                temp(objects.PixelIdxList{I(2)}) = 1;
                % Calc midline length
                cv = cHull(temp);
                [midlineTemp,DTemp] = longestLine(cv);
                % SAve results if does not exceed midline length
                if ((DTemp < maxMidlineLength) &&...
                        (length(objects.PixelIdxList{I(2)})) > minChunkSize)
                    bw = temp; midline = midlineTemp; D = DTemp;
                    % Adds third largest object if detection is too small
                    if (length(I) > 2) && (DTemp < maxMidlineLength)
                        temp = bw;
                        temp(objects.PixelIdxList{I(3)}) = 1;
                        % Calc midline length
                        cv = cHull(temp);
                        [midlineTemp,DTemp] = longestLine(cv);
                        if ((DTemp < maxMidlineLength) &&...
                                (length(objects.PixelIdxList{I(3)})) > minChunkSize)
                            bw = temp; midline = midlineTemp; D = DTemp;
                        end        
                    end
                end
            end
            
            %% Breakpoint for testing
            if i == 195
                disp('stop')
            end
            
            %% Draw Convex hull and calculate midline
            cv = cHull(bw);
            [midline,D] = longestLine(cv);
            if(midline(1,1) > midline(2,1))
                %% Make sure first popints is on left side
                midline = midline([2 1],:);
            end
            %% Get detected pixels perpendicular to midline
            n = midlineResolution; %Get numper of points to plot
            res = D/n; % interval between perpendicular lines
            % Get slope ( y1 - y2 / x1 - x2 )
            mMidline = (midline(2,2) - midline(1,2)) / (midline(2,1) - midline(1,1));
            if mMidline == 0; mMidline = 0.0000001;end;
            if mMidline == -Inf; mMidline = 999999;end;

            bMidline = midline(1,2) - mMidline*midline(1,1);
            m = -1/mMidline; % Perpendicular slope
            bStart = midline(1,2) - m*midline(1,1); %intercept of first perpendicular line
            
            % calc x/y offset for each perpendicular line segment
            dy = res/cosd(atand(mMidline)+90);
            
            % get coords of detected pixels
            [row,col,~] = find(bw);
            % Plot detected areas
            figure(1);subplot(4,1,1);hold off;imshow(bw*0.5 + bwSubtract*0.5);
            axis equal;hold on;
            plot(midline(:,1),midline(:,2),'color',[0 1 0]);
            xlim(bbox(1:2)); ylim(bbox(3:4));

            %% Initialize midline variables
            detections = [];
            transposedPoints = [];
            distanceFromMidline = [];
            
            %% Loop through each midline segment
            for j=0:n
                % plot line
                b = bStart - (j*dy);
                plot(bbox(1:2),m*bbox(1:2) + b,'color',[0 1 0])
                % Find intersecting pixels
                x = [(col - 0.5) (col + 0.5)];
                y = [(row - 0.5) (row + 0.5)];
                % predict y values for intersecting points
                yPred = m*x + b;
                % Save intersected pixels
                [idx,~,~] = find((y(:,1) <= max(yPred,[],2)) & (y(:,2) >= min(yPred,[],2)));
                %% Only execute if pixels have been detected
                if size(idx,1) > 0
                    %% Pixels which cross the adjacent line
                    detectionsTemp = [col(idx) row(idx)];
                    detectionsTemp = [mean(detectionsTemp(:,1)) mean(detectionsTemp(:,2)) ...
                        j+1 midline(1,1) midline(1,2)];
                    detections = [detections; detectionsTemp];

                    %% Snap points to adjacent lines
                    transposedPointsTemp = (detectionsTemp(:,1) + detectionsTemp(:,2)*m - m*b) / (1 + m^2);
                    transposedPointsTemp  = [transposedPointsTemp (m*transposedPointsTemp  + b)];
                    transposedPoints = [transposedPoints; transposedPointsTemp];
                    
                    %% Distance of fitted points to adjacent midline segment
                    for k = 1:(size(transposedPointsTemp,1)) 
                        % Absolute distance
                        Q1 = midline(1,:); Q2 = midline(2,:); P = transposedPointsTemp(k,:);
                        temp = abs( det([P-Q1;Q2-Q1]) )/norm(Q2-Q1);
                        % Correct negative values
                        if find(P(2) > (mMidline*P(1) + bMidline));
                            temp = temp*-1;
                        end
                        distanceFromMidline = [distanceFromMidline; temp];
                    end
                end
            end
            
            %% Add fitted spline (smooth using loess)            
            curve = smooth((detections(:,3)-1)*res,...
                distanceFromMidline,loessSmooth,'rloess');
            curve = [(detections(:,3)-1)*res curve];

            %% Interpolate missing points (cubic spline)
            temp = interp1(curve(:,1),curve(:,2),...
                (0:n)*res,'spline');                    
            curve = [((0:n)*res)' temp'];

            %% Transform nromalized points back to origional image
            xDisp= curve(:,1) * cos(atan(mMidline));
            yDisp= -curve(:,2) + (curve(:,1) * tan(atan(mMidline)));
            absXSmooth = xDisp + midline(1,1);
            absYSmooth = midline(1,2) + yDisp;
            
            %% Store results
            points = size(curve,1);
            data = [curve ...
                repmat(midline(1,1),points,1) repmat(midline(1,2),points,1 ),...
                repmat(midline(2,1),points,1 ) repmat(midline(2,2),points,1 )...
                (((1:points)-1)*res)' repmat(i,points,1) absXSmooth absYSmooth...
                repmat(mMidline,points,1)];

            %% Data columns
            % data = [fittedX(nromalized) fittedY(normalized)...
            % Midline1X Midline1Y MidlineyX Midline2Y...
            % PositionAlongMidline frame]
      
            %% Redefine rectangle area (recenter)
            % This updates the motion detection area for the following frame
            c = [((midline(1,1) + midline(2,1))/2) ((midline(1,2) + midline(2,2))/2)];
            bbox = [(c(1) - ((bbox(2) - bbox(1))/2)) ...
                (c(1) + ((bbox(2) - bbox(1))/2)) ...
                c(2) - ((bbox(4) - bbox(3))/2) ...
                c(2) + ((bbox(4) - bbox(3))/2)];

            %% Show detectons (fitted vs pixel centers)
            subplot(4,1,1);hold on;
            scatter(detections(:,1),detections(:,2),[],[0.5 0 0],'o');
            scatter(transposedPoints(:,1),transposedPoints(:,2),[],[1 0 0],'x');
            title('Midline Pixel Detection');
            %% Show normalized and fitted curve
            subplot(4,1,2); hold off;
            scatter((detections(:,3)-1)*res,distanceFromMidline,...
                [],[1 0 0],'filled','o'); 
            xlim([min(data(:,1)) max(data(:,1))]); axis equal;
            hold on; plot(data(:,1),data(:,2)); grid on;
            legend('Average Midline','Smoothed Midline');
            title('Normalized & Smoothed midline')
            %% Show image subtraction
            figure(1);subplot(4,1,3);hold off;
            imshow(greyScale(bbox(3):bbox(4),bbox(1):bbox(2)) - currentImage(bbox(3):bbox(4),bbox(1):bbox(2)))            
            title('Image subtraction (Ref image before current frame)'); hold on;
            text(10,10,...
                ['Pixels Detected: ' num2str(length(find(bw)))...
                ' Frame: ' num2str(i)...
                ' Straight Midline Length: ' num2str(D)],...
                'color',[0.8 0 0]);
            %% Show origional image overlayed with smoothing spline
            subplot(4,1,4); hold off;
            imshow(currentImage); hold on; xlim(bbox(1:2)); ylim(bbox(3:4));
            grid on;
            scatter(transposedPoints(:,1),transposedPoints(:,2),...
                [],[1 0 0],'o'); 
            plot(data(:,9),data(:,10),'color',[1 0 0],'linewidth',2);
            title('Smoothed midline overlayed on source image');
            
            if video == 1;
                 writeVideo(vw,getframe(1));
            end
            
            %% Save results
            results = [results; data];
        end
    end
    
    fw = fopen(['TrackingResuts' Individual Treatment '_' date '.csv'],'w');
    fprintf(fw,'%s \n','CurveX,CurveY,Midline1X,Midline1Y,Midline2X,Midline2Y,DistanceAlongMidline,Frame,AbsX,AbsY,MidlineSlope');
    fprintf(fw,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n',results');
    fclose(fw)
    if video == 1;
        close(vw);
    end
    disp('Analysis complete!')
 
end