    %% difference between images (return a binary image)
    function bw = backgroundSubtraction(refImage,currentImage, bbox)
        gdiff = abs(refImage - currentImage);
        thresh = 0.04;
        dims = size(currentImage);
        % Empty black image
        bw = zeros(dims(1),dims(2));
        % Insert eel bw region
        bw(bbox(3):bbox(4),bbox(1):bbox(2)) = im2bw(gdiff(bbox(3):bbox(4),bbox(1):bbox(2)),thresh);   
    end