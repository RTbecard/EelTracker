function idx = findPixelsIntersectingLine(x,y,m,b)
    % Returns the index of all pixels intersecting the line.
    % All pixels are assumed to have a width and height of 1.
    %
    % x: x-coord of pixel
    % y: y-coord of pixel
    % m: Slope of line
    % b: y-intercept of line
    
    % Define x and y ranges of pixels (width and height of 1)
    x = [(x - 0.5) (x + 0.5)];
    y = [(y - 0.5) (y + 0.5)];
    % predict y values for intersecting points
    yPred = m*x + b;
    % Save intersected pixels
    [idx,~,~] = find((y(:,1) <= max(yPred,[],2)) & (y(:,2) >= min(yPred,[],2)));
end