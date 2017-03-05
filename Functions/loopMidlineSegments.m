function [detections,transposedPoints,distanceFromMidline,mMidline,bMidline] =...
    loopMidlineSegments(P1,P2,n,D,bw)
% Calculate detected pixels along lines perpendicular to the longest axis
% of the fish
%
% P1: Start point along the major axis of the fish
% P2: End point along the major axis of the fish
% midlineResolution: Number of points to estimate along the midline
% D: Length of the longest axis of the fish

    %% Number of midline points to estimate
    res = D/n; % interval between perpendicular lines
    
    %% Slope and intercept of major axis of fish
    mMidline = (P2(2) - P1(2)) / (P2(1) - P1(1));
    if mMidline == 0; mMidline = 0.0000001;end;
    if mMidline == -Inf; mMidline = 999999;end;
    bMidline = P1(2) - mMidline*P1(1);
    
    %% Constants for adjacent lines to major axis
    m = -1/mMidline; % Perpendicular slope
    bStart = P1(2) - m*P1(1); %intercept of first perpendicular line
    % calc x/y offset for each perpendicular line segment
    dy = res/cosd(atand(mMidline)+90);

    %% Get coords of detected pixels
    [row,col,~] = find(bw);

    %% Initialize midline variables
    detections = [];
    transposedPoints = [];
    distanceFromMidline = [];

    %% Loop through each midline segment
    for j=0:n
        % plot line
        b = bStart - (j*dy);
        idx = findPixelsIntersectingLine(col,row,m,b);
        %% Only execute if pixels have been detected
        if size(idx,1) > 0
            %% Take average location of detected pixels
            detectionsTemp = [col(idx) row(idx)];
            detectionsTemp = [mean(detectionsTemp(:,1)) mean(detectionsTemp(:,2)) ...
                j+1 P1(1) P1(2)];
            detections = [detections; detectionsTemp];

            %% Snap points to adjacent lines
            [xSnap,ySnap] = snapPointsToLine(detectionsTemp(:,1),...
                detectionsTemp(:,2),m,b);
            transposedPointsTemp = [xSnap ySnap];
            transposedPoints = [transposedPoints; transposedPointsTemp];

            %% Distance of fitted points to adjacent midline segment
            for k = 1:(size(transposedPointsTemp,1)) 
                distanceFromMidline = [distanceFromMidline; ...
                    point2LineDistance(P1,P2,transposedPointsTemp(k,:))];
            end
        end
    end