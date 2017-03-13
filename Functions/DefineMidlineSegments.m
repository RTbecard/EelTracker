function [mOut,bOut] =...
    loopMidlineSegments(P1,P2,n,D)
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

    mOut = [];
    bOut = [];
    %% Loop through each midline segment
    for j=0:n
        % plot line
        bOut = [bOut; (bStart - (j*dy))];
        mOut = [mOut; m];        
    end