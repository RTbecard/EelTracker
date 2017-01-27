function Tracking(video,firstDetection)
%   startFrame = 150
%   video = 'ExampleVideo/sprint_20160722_024359_PM_20160722_025355_PM.avi';
    vr = VideoReader(video);
    
    %% Initialize results
    results = [];

    %% Load reference frame
    startFrame = 1;
    buffer = 600;  % load multiple video frames at a time (increases processing speed)
    endFrame = startFrame + buffer - 1;
    frames = read(vr,[startFrame,endFrame]);
    % The reference frame is redefined on each frame of video.  It is the
    % region around the eel n frames before or after the current frame,
    % where n = refOffset.  This is useful in case the camera moves
    % slightly during the setup, then only a portion of the video has
    % tracking errors.  This can be improved so that both the reference
    % before and after the current frame is used, and the one with the
    % least motion detection is used.
    
    %% Parameters
    % Minimum size for detections (in pixels)  Any continous detected
    % obejct with les than this many pixls is ignored
    minSize = 50;
    % How many frames ahead (or behind) the reference image will be grabbed from
    refOffset = 150;
    
    %% Function for loading video frames (current frame)
    function [frame] = loadBatch(idx)
        if idx > vr.NumberOfFrames || idx < 1
            % Return a white frame (full detection) if frame number is out
            % of bounds
           frame = currentImage;
           frame(:,:) = 255;
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

    %% Fucntion for spline
    % takes average vertical position of detections
    function[x,y] = fitspline(j,k)
        x = [];
        y = [];
        for m=unique(k)'
            idx = find(k == m);
            y = [y mean(j(idx))];
            x = [x m];
        end
    end

    %% First frame where eel is fully visible
    figure(2);
    imshow(loadBatch(firstDetection))
    disp('Select upper left and lower right region where eel is')
    [x,y] = ginput(2);
    x = int32(x);
    y = int32(y);
    close(2)
    bbox = [x(1) x(2) y(1) y(2)];
    
    %% difference between images (return a binary image)
    function bw = backgroundSubtraction(refImage)
        gdiff = abs(refImage - currentImage);
        thresh = 0.04;
        dims = size(currentImage);
        % Empty black image
        bw = zeros(dims(1),dims(2));
        % Insert eel bw region
        bw(bbox(3):bbox(4),bbox(1):bbox(2)) = im2bw(gdiff(bbox(3):bbox(4),bbox(1):bbox(2)),thresh);   
    end
    
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
        bwBefore = backgroundSubtraction(refImageBefore);
        bwAfter = backgroundSubtraction(refImageAfter);
        if length(find(bwBefore == 1)) <= length(find(bwAfter == 1)) 
            bw = bwBefore;
        else
            bw = bwAfter;
        end
        
        %% Remove small continous regions of motion
        bw = bwareaopen(bw, minSize);
        
        %% Redefine rectangle area (recenter)
        rp = regionprops(bw, 'Centroid');
        c = rp.Centroid;
        bbox = [(c(1) - ((bbox(2) - bbox(1))/2)) ...
            (c(1) + ((bbox(2) - bbox(1))/2)) ...
            c(2) - ((bbox(4) - bbox(3))/2) ...
            c(2) + ((bbox(4) - bbox(3))/2)];
        
        %% Image processing
        % remove long thin detected areas
        bw2 = bwmorph(bw, 'majority');
        
        
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