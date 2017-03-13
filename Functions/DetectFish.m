function [ref,greyScale] = DetectFish(refImageBefore,refImageAfter,...
    currentImage,bbox,thresh)

    %% Background subtract from both reference frames
	bwBefore = backgroundSubtraction(refImageBefore,currentImage, bbox,thresh);
    bwAfter = backgroundSubtraction(refImageAfter,currentImage, bbox,thresh);
    
    %% Select the subtracted image with the least detected objects
    % Assume that errors wil result in extra detected objects
    if length(find(bwBefore == 1)) <= length(find(bwAfter == 1))
        ref= bwBefore;
        greyScale = refImageBefore; 
    else
        ref = bwAfter;
        greyScale = refImageAfter; 
    end

%% difference between images (return a binary image)
    function bw = backgroundSubtraction(refImage,currentImage, bbox,thresh)
        gdiff = abs(refImage - currentImage);
        dims = size(currentImage);
        % Empty black image
        bw = zeros(dims(1),dims(2));
        % Insert eel bw region
        bw(bbox(3):bbox(4),bbox(1):bbox(2)) = ...
            im2bw(gdiff(bbox(3):bbox(4),bbox(1):bbox(2)),thresh);   
    end
end